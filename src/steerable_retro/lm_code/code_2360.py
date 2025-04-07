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
    Detects convergent synthesis strategy by checking if the route has multiple branches.
    """
    # Count the number of leaf nodes (starting materials)
    leaf_count = 0
    branch_points = 0

    def analyze_structure(node):
        nonlocal leaf_count, branch_points

        if node["type"] == "mol" and node.get("in_stock", False):
            leaf_count += 1
            return

        # Count branch points (nodes with multiple children)
        children = node.get("children", [])
        if len(children) >= 2:
            branch_points += 1

        # Recursively check children
        for child in children:
            analyze_structure(child)

    analyze_structure(route)

    # A convergent synthesis typically has multiple starting materials
    is_convergent = leaf_count >= 2
    print(
        f"Convergent synthesis analysis: {leaf_count} starting materials, {branch_points} branch points"
    )
    if is_convergent:
        print(
            f"Convergent synthesis detected with {leaf_count} starting materials and {branch_points} branch points"
        )

    return is_convergent
