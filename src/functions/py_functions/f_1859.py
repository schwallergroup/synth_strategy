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
    This function detects convergent synthesis where multiple fragments are prepared separately then coupled.
    """
    # Track branch depths
    branch_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Store the current depth
            node["depth"] = depth

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_count = len(rsmi.split(">")[0].split("."))

            # If a reaction has multiple reactants and is not at the deepest level,
            # it might be a convergent point
            if (
                reactants_count >= 2 and depth <= 2
            ):  # Focus on reactions in the first few steps
                branch_depths.append(depth)

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have multiple branches at different depths
    is_convergent = len(branch_depths) >= 1 and any(d <= 2 for d in branch_depths)

    if is_convergent:
        print(
            f"Detected convergent synthesis with branch points at depths: {branch_depths}"
        )

    return is_convergent
