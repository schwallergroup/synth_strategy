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
    # Track the number of fragments combined at each step
    fragment_combinations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count the number of distinct reactants
            num_reactants = len(reactants)
            fragment_combinations.append((depth, num_reactants))

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze the pattern of fragment combinations
    # In a linear synthesis, most steps combine only 1-2 fragments
    # In a convergent synthesis, there would be steps combining 3+ fragments

    # Sort by depth
    fragment_combinations.sort()

    # Check if most steps combine only 1-2 fragments
    small_combinations = sum(1 for _, num in fragment_combinations if num <= 2)
    large_combinations = sum(1 for _, num in fragment_combinations if num > 2)

    # If more than 80% of steps combine only 1-2 fragments, consider it linear
    is_linear = (
        small_combinations
        / (
            small_combinations + large_combinations
            if small_combinations + large_combinations > 0
            else 1
        )
    ) > 0.8

    print(f"Fragment combination pattern: {fragment_combinations}")
    print(f"Linear assembly strategy detected: {is_linear}")

    return is_linear
