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
    This function detects if the synthetic route involves pyridine-containing compounds.
    """
    pyridine_detected = False

    def dfs_traverse(node):
        nonlocal pyridine_detected

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                pyridine_pattern = Chem.MolFromSmarts("[n]1[c][c][c][c][c]1")
                if mol and mol.HasSubstructMatch(pyridine_pattern):
                    print("Pyridine structure detected")
                    pyridine_detected = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return pyridine_detected
