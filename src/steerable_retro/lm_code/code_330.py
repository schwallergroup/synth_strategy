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
    This function detects if a phenol-containing building block is used in the synthesis.
    """
    phenol_found = False

    def dfs_traverse(node):
        nonlocal phenol_found

        if node["type"] == "mol" and node.get("in_stock", False):
            # Check if this starting material contains a phenol group
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[c][OH]")):
                print("Phenol building block detected")
                phenol_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return phenol_found
