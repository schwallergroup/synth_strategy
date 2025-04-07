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
    """Helper function to determine the maximum depth of the route"""
    max_depth = [0]

    def dfs(node, depth=0):
        if depth > max_depth[0]:
            max_depth[0] = depth

        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return max_depth[0]
