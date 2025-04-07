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
    Detects if a cyclopropane ring is present throughout the synthesis,
    excluding starting materials (in_stock=True).
    """
    all_nodes_have_cyclopropane = True

    def dfs_traverse(node, depth=0):
        nonlocal all_nodes_have_cyclopropane

        # Process molecule nodes that are not starting materials
        if (
            node["type"] == "mol"
            and "smiles" in node
            and not node.get("in_stock", False)
        ):
            smiles = node["smiles"]

            # Check if this molecule has a cyclopropane ring
            has_cyclopropane = checker.check_ring("cyclopropane", smiles)

            if not has_cyclopropane:
                print(
                    f"Intermediate/product without cyclopropane found at depth {depth}: {smiles}"
                )
                all_nodes_have_cyclopropane = False

        # Continue traversal for all child nodes
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return all_nodes_have_cyclopropane
