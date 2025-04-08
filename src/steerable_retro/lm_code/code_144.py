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


def main(route, min_steps=5):
    """
    Detects if the synthesis route is a linear sequence with at least min_steps steps
    """
    step_count = 0

    def dfs_traverse(node):
        nonlocal step_count

        if node["type"] == "reaction":
            step_count += 1

        # Continue traversal - for linear synthesis, we expect only one child per reaction
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Found {step_count} linear synthesis steps")
    return step_count >= min_steps
