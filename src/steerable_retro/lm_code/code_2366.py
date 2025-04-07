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
    This function detects if the synthesis involves an isoxazole-containing fragment.
    """
    isoxazole_pattern = Chem.MolFromSmarts("c1conc1")
    isoxazole_present = False

    def dfs_traverse(node):
        nonlocal isoxazole_present

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(isoxazole_pattern):
                isoxazole_present = True
                print("Isoxazole-containing fragment detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return isoxazole_present
