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
    Detects if the synthesis follows a linear (non-convergent) strategy with at least 3 steps.
    """
    max_depth = 0
    branching_factor = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, branching_factor

        # Update max depth
        max_depth = max(max_depth, depth)

        # Count children that are reactions
        reaction_children = sum(
            1 for child in node.get("children", []) if child["type"] == "reaction"
        )

        # Update max branching factor
        branching_factor = max(branching_factor, reaction_children)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis has depth >= 3 and max branching factor <= 1
    is_linear = max_depth >= 3 and branching_factor <= 1

    if is_linear:
        print(f"Found linear synthesis strategy with depth {max_depth}")

    return is_linear
