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
    Detects if one of the key fragments in the synthesis contains a sulfonamide group.
    """
    sulfonamide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_found

        if node["type"] == "mol":
            # Check if molecule contains sulfonamide group
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])[#7]")
                if mol.HasSubstructMatch(sulfonamide_pattern):
                    print("Sulfonamide group detected in molecule at depth", depth)
                    sulfonamide_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return sulfonamide_found
