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
    This function detects a linear synthesis strategy where each step builds upon the previous.
    """
    # Track the maximum branching factor in the synthesis tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count number of non-in_stock children (reactants)
            non_stock_children = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_stock_children += 1

            # Update max branching factor
            max_branching = max(max_branching, non_stock_children)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If max_branching <= 1, it's a linear synthesis
    is_linear = max_branching <= 1
    if is_linear:
        print("Detected linear synthesis strategy")
    else:
        print(f"Detected branched synthesis with max branching factor {max_branching}")

    return is_linear
