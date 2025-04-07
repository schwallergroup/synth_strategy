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
    Detects if the route preserves a benzofuran core throughout the synthesis.
    """
    benzofuran_pattern = Chem.MolFromSmarts("[c]1[c][o][c]2[c]1[cH][cH][cH][cH]2")
    all_products_have_benzofuran = True

    def dfs_traverse(node):
        nonlocal all_products_have_benzofuran

        if node["type"] == "mol" and not node.get("in_stock", False):
            # Check if this molecule has the benzofuran core
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and not mol.HasSubstructMatch(benzofuran_pattern):
                all_products_have_benzofuran = False
                print(f"Found a product without benzofuran core: {node['smiles']}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Benzofuran core preservation: {all_products_have_benzofuran}")
    return all_products_have_benzofuran
