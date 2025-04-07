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
    This function detects if the synthetic route is a linear multistep synthesis
    with at least 4 steps.
    """
    step_count = 0

    def dfs_traverse(node):
        nonlocal step_count

        if node["type"] == "reaction":
            step_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    result = step_count >= 4
    print(f"Linear multistep synthesis strategy detected: {result} (steps: {step_count})")
    return result
