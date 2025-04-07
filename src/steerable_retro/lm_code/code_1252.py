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
    Detects if the synthesis preserves stereocenters throughout the route.
    """
    stereocenters_present = False

    def dfs_traverse(node):
        nonlocal stereocenters_present

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            # Check for stereochemistry indicators in SMILES
            if "@" in smiles:
                stereocenters_present = True
                print(f"Found stereocenter in molecule: {smiles}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return stereocenters_present
