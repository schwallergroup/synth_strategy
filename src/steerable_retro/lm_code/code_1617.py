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
    This function detects if a nitrile group is preserved throughout the synthesis.
    """
    all_intermediates_have_nitrile = True

    def dfs_traverse(node):
        nonlocal all_intermediates_have_nitrile

        if node["type"] == "mol" and node.get("smiles") and not node.get("in_stock", False):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for nitrile group
                if not mol.HasSubstructMatch(Chem.MolFromSmarts("[C]#[N]")):
                    all_intermediates_have_nitrile = False
                    print(f"Found intermediate without nitrile: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitrile group preserved throughout: {all_intermediates_have_nitrile}")
    return all_intermediates_have_nitrile
