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
    builds upon the previous without convergent steps.
    """
    # In a linear synthesis, each reaction node should have exactly one molecule child
    # that is not a starting material

    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count non-starting material molecule children
            non_starting_material_count = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_starting_material_count += 1

            # If more than one non-starting material, it's not linear
            if non_starting_material_count > 1:
                is_linear = False
                print(
                    "Found non-linear step with multiple non-starting material children"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return is_linear
