#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
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
    Detects if a specific heterocycle is preserved throughout the synthesis.
    Based on the molecules in the route, we're checking for 1,3,4-oxadiazole preservation.
    """
    # Track if the target product and all intermediate molecules have the target heterocycle
    target_product_has_ring = False
    all_intermediates_have_ring = True

    # Try different naming conventions for the oxadiazole ring
    target_rings = ["oxadiazole", "1,3,4-oxadiazole"]

    # Pattern for 1,3,4-oxadiazole in SMILES: five-membered ring with O and two N atoms
    # More specific pattern for 1,3,4-oxadiazole
    oxadiazole_pattern = Chem.MolFromSmarts("c1nnco1")

    def check_for_oxadiazole(mol_smiles):
        """Helper function to check if a molecule contains an oxadiazole ring"""
        try:
            mol = Chem.MolFromSmiles(mol_smiles)
            if mol is None:
                print(f"Failed to parse SMILES: {mol_smiles}")
                return False

            # Check if the molecule contains the target ring using checker
            for ring_name in target_rings:
                if checker.check_ring(ring_name, mol_smiles):
                    print(f"Checker found {ring_name} in: {mol_smiles}")
                    return True

            # Fallback: Use direct substructure matching if checker didn't find the ring
            if mol.HasSubstructMatch(oxadiazole_pattern):
                print(f"Pattern matching found oxadiazole in: {mol_smiles}")
                return True

            print(f"No oxadiazole found in: {mol_smiles}")
            return False

        except Exception as e:
            print(f"Error processing molecule {mol_smiles}: {str(e)}")
            # Don't fail the entire check for parsing errors
            return False

    def dfs_traverse(node, depth=0):
        nonlocal target_product_has_ring, all_intermediates_have_ring

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Skip empty SMILES
            if not mol_smiles:
                return

            # Check if this is a starting material
            is_starting_material = node.get("in_stock", False)

            # If it's the target product (depth 0)
            if depth == 0:
                target_product_has_ring = check_for_oxadiazole(mol_smiles)
                if not target_product_has_ring:
                    print(f"Target product does not contain oxadiazole: {mol_smiles}")
            # If it's an intermediate (not a starting material)
            elif not is_starting_material:
                has_ring = check_for_oxadiazole(mol_smiles)
                if not has_ring:
                    all_intermediates_have_ring = False
                    print(f"Intermediate without oxadiazole found: {mol_smiles}")
            # If it's a starting material, we don't need to check
            else:
                print(f"Skipping starting material: {mol_smiles}")

        # Recursively check all children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start the traversal from the root
    dfs_traverse(route)

    # The heterocycle is preserved if both the target product and all intermediates have it
    result = target_product_has_ring and all_intermediates_have_ring
    print(f"Target product has oxadiazole: {target_product_has_ring}")
    print(f"All intermediates have oxadiazole: {all_intermediates_have_ring}")
    print(f"Heterocycle preservation result: {result}")

    return result
