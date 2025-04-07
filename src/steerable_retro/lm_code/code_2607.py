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
    Detects if the synthesis follows a linear strategy where each reaction
    has exactly one non-commercial precursor.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count non-commercial precursors
            non_commercial_count = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_commercial_count += 1
                    # Recursively check this non-commercial precursor
                    dfs_traverse(child)

            # If more than one non-commercial precursor, it's not linear
            if non_commercial_count > 1:
                is_linear = False
                print(
                    f"Found non-linear step with {non_commercial_count} non-commercial precursors"
                )

        # For molecule nodes that aren't already processed
        elif node["type"] == "mol" and not node.get("in_stock", False):
            for child in node.get("children", []):
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if is_linear:
        print("Detected linear synthesis strategy")

    return is_linear
