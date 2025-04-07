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
    This function detects a synthetic strategy involving trifluoromethyl groups in the final product.
    """
    cf3_count = 0

    def dfs_traverse(node):
        nonlocal cf3_count

        if node["type"] == "mol" and not node.get("children", []):  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            cf3_pattern = Chem.MolFromSmarts("C(F)(F)F")
            if mol:
                cf3_count = len(mol.GetSubstructMatches(cf3_pattern))
                print(f"Found {cf3_count} trifluoromethyl groups in final product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Trifluoromethyl-containing strategy detected: {cf3_count > 0}")
    return cf3_count > 0
