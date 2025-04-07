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
    Detects a linear synthesis strategy where each step builds on a single intermediate
    (no convergent steps after initial core formation)
    """
    # Track branching factor at each step
    branching_factors = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Count number of reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            num_reactants = len([r for r in reactants if r])

            # Store branching factor (number of reactants)
            branching_factors.append((depth, num_reactants))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth
    branching_factors.sort(key=lambda x: x[0])

    # Check if linear after first step (allow first step to have multiple reactants)
    is_linear = True
    if len(branching_factors) > 1:
        for i in range(1, len(branching_factors)):
            if branching_factors[i][1] > 2:  # More than 2 reactants in later steps
                is_linear = False
                break

    print(f"Linear synthesis strategy detected: {is_linear}")
    return is_linear
