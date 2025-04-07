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
    This function detects if stereochemistry is preserved throughout the synthesis.
    """
    # Track molecules with stereochemistry at each depth
    stereo_at_depths = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            if "smiles" in node:
                try:
                    smiles = node["smiles"]
                    # Check if SMILES contains '@' which indicates stereochemistry
                    if "@" in smiles:
                        if depth not in stereo_at_depths:
                            stereo_at_depths[depth] = 0
                        stereo_at_depths[depth] += 1
                        print(f"Found stereochemistry at depth {depth}")
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If stereochemistry is present at multiple depths, it's preserved
    return len(stereo_at_depths) >= 2
