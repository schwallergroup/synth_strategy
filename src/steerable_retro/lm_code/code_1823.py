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
    Detects if the synthesis uses Boc protection strategy throughout most steps
    """
    boc_protected_intermediates = 0
    total_intermediates = 0

    def dfs_traverse(node, depth=0):
        nonlocal boc_protected_intermediates, total_intermediates

        if node["type"] == "mol" and "smiles" in node and depth > 0:  # Skip final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                total_intermediates += 1

                # Check for Boc group
                boc_pattern = Chem.MolFromSmarts("[#6]([#6])([#6])([#6])[#8][#6](=[#8])[#7]")
                if mol.HasSubstructMatch(boc_pattern):
                    boc_protected_intermediates += 1
                    print(f"Found Boc-protected intermediate at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if more than half of intermediates are Boc-protected
    if total_intermediates > 0:
        ratio = boc_protected_intermediates / total_intermediates
        print(
            f"Boc protection ratio: {ratio} ({boc_protected_intermediates}/{total_intermediates})"
        )
        return ratio >= 0.5

    return False
