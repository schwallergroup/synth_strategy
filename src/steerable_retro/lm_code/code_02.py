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
    Detects a synthetic strategy involving compounds containing dimethylamine functional group.
    """
    has_dimethylamine = False

    def dfs_traverse(node):
        nonlocal has_dimethylamine

        if node["type"] == "mol":
            # Check for dimethylamine group in the molecule
            smiles = node.get("smiles", "")
            if "[N]([C][H][H][H])([C][H][H][H])" in smiles or "N(C)C" in smiles:
                has_dimethylamine = True
                print(f"Detected dimethylamine group: {smiles}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Contains dimethylamine group: {has_dimethylamine}")

    return has_dimethylamine
