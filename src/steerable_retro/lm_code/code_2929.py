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
    This function detects a synthetic strategy where a thiophene core is preserved
    throughout the synthesis. The thiophene ring should be present in the final product
    and all intermediate products, but not necessarily in reagents or starting materials.
    """
    # Track if the main synthetic pathway maintains a thiophene core
    thiophene_in_main_pathway = True

    def dfs_traverse(node, depth=0, is_main_product=True):
        nonlocal thiophene_in_main_pathway

        # Only check molecules that are in the main product line
        if node["type"] == "mol" and node.get("smiles") and is_main_product:
            mol_smiles = node["smiles"]
            is_starting_material = node.get("in_stock", False)

            # Check if molecule contains thiophene, but only if it's not a starting material
            if not is_starting_material:
                try:
                    has_thiophene = checker.check_ring("thiophene", mol_smiles)
                    if not has_thiophene:
                        print(
                            f"Intermediate at depth {depth} without thiophene found: {mol_smiles}"
                        )
                        thiophene_in_main_pathway = False
                    else:
                        print(f"Found thiophene in molecule at depth {depth}: {mol_smiles}")
                except Exception as e:
                    print(f"Error checking thiophene in molecule: {mol_smiles}, Error: {e}")
            else:
                print(f"Skipping starting material check at depth {depth}: {mol_smiles}")

        # Process children nodes
        if node["type"] == "reaction" and node.get("children"):
            # In retrosynthetic traversal, the first child is the product
            for i, child in enumerate(node.get("children", [])):
                # Only the first child is considered the main product in retrosynthetic direction
                dfs_traverse(child, depth + 1, is_main_product=(i == 0))
        else:
            # Process children for molecule nodes (maintaining the is_main_product flag)
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1, is_main_product)

    # Start traversal with the target molecule
    try:
        # Check if the target molecule contains thiophene
        if route["type"] == "mol" and route.get("smiles"):
            target_has_thiophene = checker.check_ring("thiophene", route["smiles"])
            if not target_has_thiophene:
                print(f"Target molecule does not contain thiophene: {route['smiles']}")
                return False
            else:
                print(f"Target molecule contains thiophene: {route['smiles']}")

        # Traverse the synthesis route
        dfs_traverse(route)

        return thiophene_in_main_pathway
    except Exception as e:
        print(f"Error in thiophene_based_synthesis_strategy: {e}")
        return False
