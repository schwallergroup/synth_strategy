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
    Detects if Boc protection is maintained throughout the synthesis.
    """
    boc_counts_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                boc_pattern = Chem.MolFromSmarts("[N][C$(C=O)][O][C]([C])([C])[C]")
                if mol.HasSubstructMatch(boc_pattern):
                    boc_counts_by_depth[depth] = boc_counts_by_depth.get(depth, 0) + 1
                    print(f"Found Boc protection at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if Boc protection is present at multiple depths
    return len(boc_counts_by_depth) >= 2
