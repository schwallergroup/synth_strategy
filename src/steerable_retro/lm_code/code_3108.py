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
    Detects if the synthesis involves a pyrazole heterocycle.
    """
    has_pyrazole = False

    def dfs_traverse(node):
        nonlocal has_pyrazole

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if not mol:
                return

            # Check for pyrazole pattern
            pyrazole_patt = Chem.MolFromSmarts("[n]1[n][c,n][c,n][c,n]1")
            if mol.HasSubstructMatch(pyrazole_patt):
                has_pyrazole = True
                print("Detected pyrazole heterocycle")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_pyrazole
