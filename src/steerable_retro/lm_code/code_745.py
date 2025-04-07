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
    This function detects if the synthesis begins with an ester derivative.
    """
    has_ester_start = False

    def dfs_traverse(node, depth=0):
        nonlocal has_ester_start

        if node["type"] == "mol" and node.get("in_stock", False):
            ester_pattern = Chem.MolFromSmarts("[#6][C;$(C=O)][O][#6]")
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol and mol.HasSubstructMatch(ester_pattern):
                has_ester_start = True
                print(f"Ester starting material detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_ester_start
