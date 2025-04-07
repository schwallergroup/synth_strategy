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
    This function detects a synthetic strategy that utilizes an aromatic amine
    as a key precursor.
    """
    aromatic_amine_found = False

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_amine_found

        if node["type"] == "mol" and depth > 1:  # Look for precursors (higher depth)
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for aromatic amine
                aromatic_amine_pattern = Chem.MolFromSmarts("[c][NH2]")
                if mol.HasSubstructMatch(aromatic_amine_pattern):
                    aromatic_amine_found = True
                    print(f"Aromatic amine precursor found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if aromatic_amine_found:
        print("Aromatic amine precursor strategy detected")

    return aromatic_amine_found
