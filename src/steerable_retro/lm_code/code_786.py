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
    Detects formation of a fused thieno-pyridine heterocyclic system.
    """
    thieno_pyridine_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal thieno_pyridine_formation_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains thieno-pyridine system
            product_mol = Chem.MolFromSmiles(product)

            # Check for thieno-pyridine in product using the checker function
            has_thieno_pyridine = checker.check_ring("thiophene", product) and checker.check_ring(
                "pyridine", product
            )

            if has_thieno_pyridine and product_mol:
                print(f"Product contains thiophene and pyridine rings: {product}")

                # Verify that the thiophene and pyridine rings are fused
                # by checking if they share atoms in the product
                thiophene_indices = checker.get_ring_atom_indices("thiophene", product)
                pyridine_indices = checker.get_ring_atom_indices("pyridine", product)

                if thiophene_indices and pyridine_indices:
                    # Check if any atoms are shared between the rings (fusion)
                    fused_system = False
                    for thiophene_atoms in thiophene_indices:
                        for pyridine_atoms in pyridine_indices:
                            # Convert to sets for intersection check
                            thiophene_set = set(thiophene_atoms)
                            pyridine_set = set(pyridine_atoms)

                            if thiophene_set.intersection(pyridine_set):
                                fused_system = True
                                print(f"Found fused thieno-pyridine system in product: {product}")
                                break
                        if fused_system:
                            break

                    if fused_system:
                        # Check if reactants don't have the complete fused system
                        has_fused_system_in_reactants = False

                        for reactant in reactants:
                            reactant_mol = Chem.MolFromSmiles(reactant)
                            if reactant_mol:
                                # Check if reactant has both rings and they're fused
                                if checker.check_ring("thiophene", reactant) and checker.check_ring(
                                    "pyridine", reactant
                                ):

                                    # Check if they're fused in the reactant
                                    r_thiophene_indices = checker.get_ring_atom_indices(
                                        "thiophene", reactant
                                    )
                                    r_pyridine_indices = checker.get_ring_atom_indices(
                                        "pyridine", reactant
                                    )

                                    if r_thiophene_indices and r_pyridine_indices:
                                        for r_thiophene_atoms in r_thiophene_indices:
                                            for r_pyridine_atoms in r_pyridine_indices:
                                                r_thiophene_set = set(r_thiophene_atoms)
                                                r_pyridine_set = set(r_pyridine_atoms)

                                                if r_thiophene_set.intersection(r_pyridine_set):
                                                    has_fused_system_in_reactants = True
                                                    print(
                                                        f"Fused system already exists in reactant: {reactant}"
                                                    )
                                                    break
                                            if has_fused_system_in_reactants:
                                                break

                        if not has_fused_system_in_reactants:
                            # Check if this is a ring-forming reaction
                            is_ring_forming = False

                            # Check for common ring-forming reaction types
                            if (
                                checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                                or checker.check_reaction("{benzothiophene}", rsmi)
                                or checker.check_reaction("{thiazole}", rsmi)
                                or checker.check_reaction("{indole}", rsmi)
                                or checker.check_reaction(
                                    "{benzimidazole_derivatives_aldehyde}", rsmi
                                )
                                or checker.check_reaction("{benzoxazole_arom-aldehyde}", rsmi)
                            ):
                                is_ring_forming = True
                                print(f"Detected ring-forming reaction: {rsmi}")

                            # Check if reactants contain components that could form thieno-pyridine
                            has_thiophene = False
                            has_pyridine = False
                            has_thiophene_precursor = False
                            has_pyridine_precursor = False

                            for reactant in reactants:
                                if checker.check_ring("thiophene", reactant):
                                    has_thiophene = True
                                    print(f"Reactant contains thiophene: {reactant}")
                                if checker.check_ring("pyridine", reactant):
                                    has_pyridine = True
                                    print(f"Reactant contains pyridine: {reactant}")
                                # Check for precursors that could form these rings
                                if (
                                    checker.check_fg("Thioamide", reactant)
                                    or checker.check_fg("Aliphatic thiol", reactant)
                                    or checker.check_fg("Aromatic thiol", reactant)
                                ):
                                    has_thiophene_precursor = True
                                    print(f"Reactant contains thiophene precursor: {reactant}")
                                if (
                                    checker.check_fg("Nitrile", reactant)
                                    or checker.check_fg("Primary amine", reactant)
                                    or checker.check_fg("Secondary amine", reactant)
                                ):
                                    has_pyridine_precursor = True
                                    print(f"Reactant contains pyridine precursor: {reactant}")

                            # If we have separate components that become fused, it's a formation
                            if (
                                is_ring_forming
                                or (has_thiophene and has_pyridine)
                                or (has_thiophene and has_pyridine_precursor)
                                or (has_thiophene_precursor and has_pyridine)
                                or (has_thiophene_precursor and has_pyridine_precursor)
                            ):
                                print(f"Thieno-pyridine formation detected in reaction: {rsmi}")
                                thieno_pyridine_formation_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return thieno_pyridine_formation_detected
