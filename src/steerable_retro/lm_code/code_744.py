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
    This function detects if the synthesis maintains fluorinated aromatic rings.
    """
    fluoro_aromatic_counts = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node.get("smiles"):
            fluoro_aromatic_pattern = Chem.MolFromSmarts("[c][F]")
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol and mol.HasSubstructMatch(fluoro_aromatic_pattern):
                fluoro_aromatic_counts.append(depth)
                print(f"Fluorinated aromatic detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if fluorinated aromatics are present at multiple depths (maintained throughout)
    if len(fluoro_aromatic_counts) >= 2:
        print("Fluorinated aromatics maintained throughout synthesis")
        return True

    return False
