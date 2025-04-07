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
    This function detects if the synthesis uses an alkyne as a rigid linker
    between aromatic systems in the final product.
    """
    has_alkyne_linker = False

    def dfs_traverse(node, depth=0):
        nonlocal has_alkyne_linker

        if node["type"] == "mol" and depth == 0:  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for alkyne connecting two aromatic systems
                alkyne_linker_pattern = Chem.MolFromSmarts("[c,n][C]#[C][c,n]")
                if mol.HasSubstructMatch(alkyne_linker_pattern):
                    has_alkyne_linker = True
                    print("Alkyne linker detected in final product")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_alkyne_linker
