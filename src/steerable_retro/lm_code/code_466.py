#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    This function detects a comprehensive synthetic strategy involving:
    1. Formation of sulfonamide-containing aromatic heterocycle
    2. Reduction of aromatic N-heterocycle to saturated heterocycle
    3. Late-stage N-alkylation to join complex fragments
    """
    # Track each component of the strategy
    has_sulfonamide_formation = False
    has_heterocycle_reduction = False
    has_late_n_alkylation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_sulfonamide_formation, has_heterocycle_reduction, has_late_n_alkylation

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonamide formation
                if checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                ):
                    print(f"Detected sulfonamide formation at depth {depth}")
                    has_sulfonamide_formation = True
                elif any(
                    checker.check_fg("Sulfonyl halide", r) for r in reactants
                ) and checker.check_fg("Sulfonamide", product):
                    print(f"Detected sulfonamide formation from sulfonyl halide at depth {depth}")
                    has_sulfonamide_formation = True

                # Check for heterocycle reduction
                # Look for aromatic heterocycles in reactants and saturated heterocycles in product
                aromatic_heterocycles = [
                    "pyridine",
                    "pyrrole",
                    "imidazole",
                    "pyrazole",
                    "oxazole",
                    "thiazole",
                ]
                saturated_heterocycles = [
                    "piperidine",
                    "pyrrolidine",
                    "piperazine",
                    "morpholine",
                    "thiomorpholine",
                ]

                for reactant in reactants:
                    for arom_het in aromatic_heterocycles:
                        if checker.check_ring(arom_het, reactant):
                            for sat_het in saturated_heterocycles:
                                if checker.check_ring(sat_het, product):
                                    print(
                                        f"Detected heterocycle reduction: {arom_het} to {sat_het} at depth {depth}"
                                    )
                                    has_heterocycle_reduction = True

                # Check for late-stage N-alkylation (depth 0 or 1)
                if depth <= 1:
                    # Check for N-alkylation reactions
                    if (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("Alkylation of amines", rsmi)
                    ):

                        # Check if at least one reactant is complex (has many atoms)
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                        reactant_heavy_atoms = [
                            sum(1 for atom in r.GetAtoms() if atom.GetAtomicNum() > 1)
                            for r in reactant_mols
                            if r
                        ]

                        if any(count >= 8 for count in reactant_heavy_atoms):
                            print(
                                f"Detected late-stage N-alkylation joining complex fragments at depth {depth}"
                            )
                            has_late_n_alkylation = True
                    # Fallback check for N-alkylation if reaction type check fails
                    elif not has_late_n_alkylation:
                        reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                        product_mol = Chem.MolFromSmiles(product) if product else None

                        if product_mol and reactant_mols and len(reactant_mols) >= 2:
                            has_amine = any(
                                checker.check_fg("Secondary amine", r)
                                or checker.check_fg("Primary amine", r)
                                for r in reactants
                            )
                            has_alkyl_halide = any(
                                checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Tertiary halide", r)
                                for r in reactants
                            )
                            has_tertiary_amine = checker.check_fg("Tertiary amine", product)

                            if has_amine and has_alkyl_halide and has_tertiary_amine:
                                reactant_heavy_atoms = [
                                    sum(1 for atom in r.GetAtoms() if atom.GetAtomicNum() > 1)
                                    for r in reactant_mols
                                    if r
                                ]
                                if any(count >= 8 for count in reactant_heavy_atoms):
                                    print(
                                        f"Detected late-stage N-alkylation via FG check at depth {depth}"
                                    )
                                    has_late_n_alkylation = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Print final detection status
    print(f"Sulfonamide formation: {has_sulfonamide_formation}")
    print(f"Heterocycle reduction: {has_heterocycle_reduction}")
    print(f"Late N-alkylation: {has_late_n_alkylation}")

    # The strategy is present if all three components are detected
    return has_sulfonamide_formation and has_heterocycle_reduction and has_late_n_alkylation
