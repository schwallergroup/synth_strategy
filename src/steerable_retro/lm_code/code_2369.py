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
    This function detects if the synthesis follows a linear rather than convergent approach.
    """
    max_branching = 0

    def count_children(node):
        if node["type"] == "reaction":
            return len(node.get("children", []))
        return 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            branching = count_children(node)
            max_branching = max(max_branching, branching)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If max branching is <= 2, it's likely a linear synthesis
    is_linear = max_branching <= 2
    print(f"Linear synthesis strategy: {is_linear} (max branching: {max_branching})")
    return is_linear
