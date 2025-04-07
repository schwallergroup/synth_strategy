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
    This function detects if the synthesis follows a linear fragment assembly strategy
    rather than a convergent approach.
    """
    reaction_count = 0
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, max_depth

        if node["type"] == "reaction":
            reaction_count += 1
            max_depth = max(max_depth, depth)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If max_depth is close to reaction_count, it suggests a linear synthesis
    is_linear = max_depth + 1 >= reaction_count
    if is_linear:
        print(
            f"Linear synthesis detected with {reaction_count} reactions and max depth {max_depth}"
        )

    return is_linear
