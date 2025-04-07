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


def main(route, threshold=3):
    """
    Detects if the final product in the synthesis route has multiple methoxy groups
    """
    methoxy_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal methoxy_count

        if node["type"] == "mol" and depth == 0:  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                methoxy_patt = Chem.MolFromSmarts("[O][CH3]")
                methoxy_count = len(mol.GetSubstructMatches(methoxy_patt))
                print(f"Found {methoxy_count} methoxy groups in final product")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return methoxy_count >= threshold
