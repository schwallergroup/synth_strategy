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
    Checks if the synthetic route is predominantly linear.
    A linear synthesis typically has a single branch at each step.
    """
    # Count total reactions and branching points
    total_reactions = 0
    branching_points = 0

    def dfs_traverse(node):
        nonlocal total_reactions, branching_points

        if node["type"] == "reaction":
            total_reactions += 1
            # Count reactants (children) for this reaction
            reactant_count = sum(
                1 for child in node.get("children", []) if child["type"] == "mol"
            )
            if reactant_count > 1:
                branching_points += 1

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Calculate branching ratio - lower is more linear
    if total_reactions > 0:
        branching_ratio = branching_points / total_reactions
        is_linear = branching_ratio <= 0.6  # Allow up to 60% branching
        print(
            f"Linearity check: {branching_points} branching points out of {total_reactions} reactions (ratio: {branching_ratio:.2f})"
        )
        return is_linear
    else:
        print("No reactions found in route")
        return False
