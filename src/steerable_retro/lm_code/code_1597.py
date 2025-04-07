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
    Detects if the final product contains an isoxazole heterocycle.
    """
    has_isoxazole = False

    def dfs_traverse(node, depth=0):
        nonlocal has_isoxazole

        if node["type"] == "mol" and depth == 0:  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if not mol:
                return

            # Check for isoxazole pattern
            isoxazole_pattern = Chem.MolFromSmarts("c1oncc1")
            if mol.HasSubstructMatch(isoxazole_pattern):
                has_isoxazole = True
                print("Detected isoxazole in final product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Isoxazole detected: {has_isoxazole}")
    return has_isoxazole
