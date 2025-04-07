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
    This function detects if the synthesis maintains a pyridine scaffold throughout.
    """
    all_intermediates_have_pyridine = True

    def dfs_traverse(node):
        nonlocal all_intermediates_have_pyridine

        if node["type"] == "mol" and "smiles" in node and not node.get("in_stock", False):
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            # Pyridine pattern (allowing for N-oxide)
            pyridine_pattern = Chem.MolFromSmarts("[n;r6]")

            if mol and not mol.HasSubstructMatch(pyridine_pattern):
                all_intermediates_have_pyridine = False
                print(f"Intermediate without pyridine scaffold detected: {smiles}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return all_intermediates_have_pyridine
