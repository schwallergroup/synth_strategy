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
    This function detects if the synthesis involves a mesityl (2,4,6-trimethylphenyl) group.
    """
    mesityl_found = False

    def dfs_traverse(node):
        nonlocal mesityl_found

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                # Mesityl group SMARTS pattern
                mesityl_pattern = Chem.MolFromSmarts("[c]1[c]([CH3])[c][c]([CH3])[c][c]1[CH3]")
                if mol and mol.HasSubstructMatch(mesityl_pattern):
                    mesityl_found = True
                    print("Mesityl group found in molecule")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return mesityl_found
