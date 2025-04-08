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


def main(route):
    """
    This function detects a linear synthesis strategy with sequential transformations
    rather than convergent fragment coupling.
    """
    # Count the number of reactions and check for branching
    reaction_count = 0
    max_children_per_node = 0

    def dfs_traverse(node):
        nonlocal reaction_count, max_children_per_node

        if node["type"] == "reaction":
            reaction_count += 1
            children_count = len(node.get("children", []))
            max_children_per_node = max(max_children_per_node, children_count)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # A linear synthesis typically has:
    # 1. Multiple reactions (>3)
    # 2. No more than 2 children per node (typically just 1, but allowing 2 for cases where a reagent is also represented)
    is_linear = reaction_count >= 3 and max_children_per_node <= 2

    if is_linear:
        print(f"Detected linear synthesis with {reaction_count} sequential steps")

    return is_linear
