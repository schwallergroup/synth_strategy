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
    Detects if the synthesis follows a linear pattern where each reaction
    builds on the product of the previous reaction.
    """
    # Count total reactions and branching points
    total_reactions = 0
    branching_points = 0

    def dfs_traverse(node, depth=0):
        nonlocal total_reactions, branching_points

        if node["type"] == "reaction":
            total_reactions += 1

            # Count children that are reactions (branching)
            reaction_children = sum(
                1 for child in node.get("children", []) if child["type"] == "reaction"
            )
            if reaction_children > 1:
                branching_points += 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # If there are no branching points and multiple reactions, it's a linear synthesis
    is_linear = (total_reactions > 1) and (branching_points == 0)
    if is_linear:
        print(f"Found linear synthesis with {total_reactions} reactions")

    return is_linear
