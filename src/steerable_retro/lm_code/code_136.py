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
    This function detects if the synthesis maintains stereochemistry
    throughout the route.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            # Check for stereochemistry indicators in SMILES
            if "@" in smiles:
                print(f"Found stereochemistry at depth {depth}: {smiles}")
                stereocenters_by_depth[depth] = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if stereocenters are present at multiple depths
    return len(stereocenters_by_depth) >= 2
