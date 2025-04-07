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
    This function detects if the synthesis route maintains a benzodioxane scaffold
    throughout the synthesis.
    """
    benzodioxane_present = False

    def dfs_traverse(node):
        nonlocal benzodioxane_present

        if node["type"] == "mol" and "smiles" in node:
            # Check for benzodioxane scaffold
            benzodioxane_pattern = Chem.MolFromSmarts("c1ccc2c(c1)OCC[C]2")

            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(benzodioxane_pattern):
                print("Benzodioxane scaffold detected")
                benzodioxane_present = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return benzodioxane_present
