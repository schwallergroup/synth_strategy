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
    This function detects a synthetic strategy involving late-stage lactam formation
    to create a heterocyclic ring system.
    """
    # Track if we found the strategy
    found_strategy = False
    # Define what "late-stage" means (low depth in retrosynthesis tree)
    LATE_STAGE_THRESHOLD = 2

    def dfs_traverse(node, depth=0):
        nonlocal found_strategy

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (low depth)
            if depth <= LATE_STAGE_THRESHOLD:
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Convert to RDKit molecules
                try:
                    reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                    product = Chem.MolFromSmiles(product_smiles)

                    if product and all(reactants):
                        # Count rings in reactants and product
                        reactant_rings = sum([len(r.GetRingInfo().AtomRings()) for r in reactants])
                        product_rings = len(product.GetRingInfo().AtomRings())

                        # Check if this is a ring formation reaction
                        if product_rings > reactant_rings:
                            print(f"Ring formation detected: {reactant_rings} → {product_rings}")

                            # Check for amide/lactam in product - include tertiary amides and specific lactam rings
                            has_amide_in_product = (
                                checker.check_fg("Secondary amide", product_smiles)
                                or checker.check_fg("Primary amide", product_smiles)
                                or checker.check_fg("Tertiary amide", product_smiles)
                                or checker.check_ring("pyrrolidone", product_smiles)
                            )

                            # Additional check for lactam structures in rings
                            if not has_amide_in_product:
                                # Look for C=O and N in the same ring
                                product_mol = Chem.MolFromSmiles(product_smiles)
                                if product_mol:
                                    ring_info = product_mol.GetRingInfo()
                                    for ring in ring_info.AtomRings():
                                        ring_atoms = set(ring)
                                        # Find carbonyl carbons in the ring
                                        carbonyl_carbons = []
                                        for atom_idx in ring_atoms:
                                            atom = product_mol.GetAtomWithIdx(atom_idx)
                                            if atom.GetSymbol() == "C":
                                                for bond in atom.GetBonds():
                                                    other_atom = bond.GetOtherAtom(atom)
                                                    if (
                                                        other_atom.GetSymbol() == "O"
                                                        and bond.GetBondType()
                                                        == Chem.BondType.DOUBLE
                                                    ):
                                                        carbonyl_carbons.append(atom_idx)
                                                        break

                                        # Check if any carbonyl carbon is adjacent to a nitrogen in the ring
                                        for c_idx in carbonyl_carbons:
                                            c_atom = product_mol.GetAtomWithIdx(c_idx)
                                            for bond in c_atom.GetBonds():
                                                other_atom = bond.GetOtherAtom(c_atom)
                                                if (
                                                    other_atom.GetSymbol() == "N"
                                                    and other_atom.GetIdx() in ring_atoms
                                                ):
                                                    has_amide_in_product = True
                                                    print(f"Detected lactam structure in ring")
                                                    break
                                            if has_amide_in_product:
                                                break

                            print(f"Has amide/lactam in product: {has_amide_in_product}")

                            # Check for amide formation reactions
                            is_amide_formation = (
                                checker.check_reaction(
                                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                                )
                                or checker.check_reaction(
                                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                                    rsmi,
                                )
                                or checker.check_reaction(
                                    "Intramolecular amination (heterocycle formation)", rsmi
                                )
                                or checker.check_reaction(
                                    "Carboxylic acid with primary amine to amide", rsmi
                                )
                                or checker.check_reaction(
                                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                                    rsmi,
                                )
                                or checker.check_reaction(
                                    "Intramolecular transesterification/Lactone formation", rsmi
                                )
                                or checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                                or checker.check_reaction("Ester with primary amine to amide", rsmi)
                                or checker.check_reaction(
                                    "Ester with secondary amine to amide", rsmi
                                )
                            )
                            print(f"Is amide formation reaction: {is_amide_formation}")

                            # Check for cyclization reactions that might form lactams
                            is_cyclization = (
                                checker.check_reaction(
                                    "Intramolecular amination (heterocycle formation)", rsmi
                                )
                                or checker.check_reaction("Formation of NOS Heterocycles", rsmi)
                                or checker.check_reaction(
                                    "Intramolecular transesterification/Lactone formation", rsmi
                                )
                            )

                            # Check for reactant functional groups that could form lactams
                            has_potential_lactam_reactants = False
                            for reactant_smile in reactants_smiles:
                                # Check for carbonyl groups
                                has_carbonyl = (
                                    checker.check_fg("Carboxylic acid", reactant_smile)
                                    or checker.check_fg("Ester", reactant_smile)
                                    or checker.check_fg("Acyl halide", reactant_smile)
                                    or checker.check_fg("Amide", reactant_smile)
                                )

                                # Check for nitrogen nucleophiles
                                has_nitrogen = (
                                    checker.check_fg("Primary amine", reactant_smile)
                                    or checker.check_fg("Secondary amine", reactant_smile)
                                    or checker.check_fg("Primary amide", reactant_smile)
                                    or checker.check_fg("Secondary amide", reactant_smile)
                                )

                                if has_carbonyl and has_nitrogen:
                                    has_potential_lactam_reactants = True
                                    print(
                                        f"Found potential lactam-forming groups in reactant: {reactant_smile}"
                                    )
                                    break

                            # Determine if this is a lactam cyclization
                            is_lactam_cyclization = (
                                product_rings > reactant_rings
                                and has_amide_in_product
                                and (
                                    is_amide_formation
                                    or is_cyclization
                                    or has_potential_lactam_reactants
                                )
                            )

                            if is_lactam_cyclization:
                                print(f"Found late-stage lactam cyclization at depth {depth}")
                                found_strategy = True
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_strategy
