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
    This function detects if trifluoromethyl-containing building blocks are used in the synthesis.
    """
    cf3_found = False

    def dfs_traverse(node):
        nonlocal cf3_found

        if node["type"] == "mol" and node.get("in_stock", False):
            # Check if this starting material contains a CF3 group
            smiles = node.get("smiles", "")
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                cf3_pattern = Chem.MolFromSmarts("[#6]-[#6]([F])([F])[F]")
                if mol.HasSubstructMatch(cf3_pattern):
                    cf3_found = True
                    print(f"Found CF3-containing building block: {smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return cf3_found
