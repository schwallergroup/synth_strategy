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
    Detects if one of the key fragments in the synthesis contains a dimethoxy-substituted aromatic system.
    """
    dimethoxy_found = False

    def dfs_traverse(node, depth=0):
        nonlocal dimethoxy_found

        if node["type"] == "mol":
            # Check if molecule contains dimethoxy aromatic pattern
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                dimethoxy_pattern = Chem.MolFromSmarts(
                    "[#6]1:[#6]:[#6]([#8][#6]):[#6]:[#6]([#8][#6]):[#6]:1"
                )
                if mol.HasSubstructMatch(dimethoxy_pattern):
                    print("Dimethoxy aromatic pattern detected in molecule at depth", depth)
                    dimethoxy_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return dimethoxy_found
