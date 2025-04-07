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
    This function detects if a benzyl group is maintained throughout the synthesis.
    """
    benzyl_pattern = Chem.MolFromSmarts("[#6]-c1ccccc1")
    benzyl_present = False

    def dfs_traverse(node):
        nonlocal benzyl_present

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(benzyl_pattern):
                benzyl_present = True
                print(f"Benzyl group found in molecule: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return benzyl_present
