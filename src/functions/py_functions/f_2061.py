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
    This function detects if a synthetic route follows a linear synthesis strategy
    without convergent steps.

    A linear synthesis has no convergent steps, meaning at each reaction step,
    there is at most one non-in_stock reactant that continues the synthesis path.
    """
    # Track the maximum branching factor in the synthesis tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction" and "children" in node:
            # Count number of non-in_stock molecule children
            # These represent the main synthetic pathway, not reagents/catalysts
            non_stock_reactants = 0
            for child in node["children"]:
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_stock_reactants += 1

            max_branching = max(max_branching, non_stock_reactants)
            print(f"Reaction node: non-stock reactants = {non_stock_reactants}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # If max_branching <= 1, it's a linear synthesis
    # If max_branching > 1, it has convergent steps
    is_linear = max_branching <= 1
    print(
        f"Linear synthesis strategy detected: {is_linear} (max branching: {max_branching})"
    )
    return is_linear
