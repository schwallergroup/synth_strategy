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
    Detects if the synthesis involves assembling multiple fragments (3+).
    Counts unique fragments that are combined throughout the synthesis.
    """
    fragments = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node.get("in_stock", False):
            # This is a starting material
            fragments.add(node["smiles"])

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have 3 or more unique fragments
    if len(fragments) >= 3:
        print(f"Found multi-fragment assembly with {len(fragments)} fragments")
        return True
    return False
