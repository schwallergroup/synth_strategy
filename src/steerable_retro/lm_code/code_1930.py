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
    This function detects if the synthesis follows a linear strategy where each step
    builds on the previous one, rather than a convergent approach.
    """
    # In a linear synthesis, most reactions have only one non-commercial reactant
    # and the tree structure is mostly linear

    reaction_count = 0
    linear_reaction_count = 0
    max_branching = 0

    def count_non_commercial_reactants(node):
        count = 0
        for child in node.get("children", []):
            if child["type"] == "mol" and not child.get("in_stock", False):
                count += 1
        return count

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, linear_reaction_count, max_branching

        if node["type"] == "reaction":
            reaction_count += 1

            # Count non-commercial reactants
            non_commercial_count = count_non_commercial_reactants(node)

            # Update max branching
            max_branching = max(max_branching, non_commercial_count)

            # If only one non-commercial reactant, it's a linear step
            if non_commercial_count <= 1:
                linear_reaction_count += 1
                print(f"Detected linear reaction step at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Calculate linearity ratio
    linearity_ratio = linear_reaction_count / reaction_count if reaction_count > 0 else 0

    print(f"Linearity ratio: {linearity_ratio}, Max branching: {max_branching}")

    # Return True if the synthesis is predominantly linear
    return linearity_ratio >= 0.7 and max_branching <= 1
