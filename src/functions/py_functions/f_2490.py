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
    Detects if the final product contains both sulfonamide and thiazole groups.
    """
    contains_both_groups = False

    def dfs_traverse(node, depth=0):
        nonlocal contains_both_groups

        # Check only the final product (depth 0)
        if depth == 0 and node["type"] == "mol":
            print(f"Checking final product: {node['smiles']}")

            # Use checker functions to detect functional groups
            has_sulfonamide = checker.check_fg("Sulfonamide", node["smiles"])
            has_thiazole = checker.check_ring("thiazole", node["smiles"])

            print(f"Has sulfonamide: {has_sulfonamide}")
            print(f"Has thiazole: {has_thiazole}")

            if has_sulfonamide and has_thiazole:
                contains_both_groups = True
                print("Final product contains both sulfonamide and thiazole groups")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return contains_both_groups
