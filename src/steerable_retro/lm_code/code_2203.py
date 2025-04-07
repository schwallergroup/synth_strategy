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
    This function detects if the synthesis route is convergent (multiple fragments combined).
    """
    # Track the number of major branches in the synthesis tree
    branch_count = 0

    def count_branches(node):
        if node["type"] == "mol" and node.get("in_stock", False):
            # This is a starting material
            return 1

        if not node.get("children", []):
            return 0

        # For reaction nodes, count branches from children
        total_branches = 0
        for child in node.get("children", []):
            total_branches += count_branches(child)

        return total_branches

    # Count branches from the root
    branch_count = count_branches(route)

    # A convergent synthesis typically has 3 or more branches
    is_convergent = branch_count >= 3
    if is_convergent:
        print(f"Convergent synthesis detected with {branch_count} branches")

    return is_convergent
