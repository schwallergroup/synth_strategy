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
    This function detects the use of an alkyne as a rigid linker between aromatic rings.
    """
    has_alkyne_linker = False

    def dfs_traverse(node):
        nonlocal has_alkyne_linker

        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for aromatic-alkyne-aromatic pattern
                alkyne_linker_pattern = Chem.MolFromSmarts("a-C#C-a")
                if mol.HasSubstructMatch(alkyne_linker_pattern):
                    print(f"Detected alkyne linker between aromatic rings: {node['smiles']}")
                    has_alkyne_linker = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    return has_alkyne_linker
