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
    This function detects a synthetic strategy involving a benzoxazinone fragment
    in the final product.
    """
    benzoxazinone_found = False

    def dfs_traverse(node, depth=0):
        nonlocal benzoxazinone_found

        if node["type"] == "mol" and depth == 0:  # Final product
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                benzoxazinone_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6]2[#6]([#6][#6][#6]1)[#8][#6][#6](=[O])[#7]2"
                )

                if mol and benzoxazinone_pattern and mol.HasSubstructMatch(benzoxazinone_pattern):
                    print("Found benzoxazinone fragment in final product")
                    benzoxazinone_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return benzoxazinone_found
