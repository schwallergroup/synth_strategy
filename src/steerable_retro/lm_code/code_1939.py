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
    This function detects if a diaryl ether motif is preserved throughout the synthesis.
    It checks if the diaryl ether motif is preserved in the target molecule and all intermediates
    after it first appears in the synthesis route.
    """
    # Track if diaryl ether appears and is then preserved
    motif_appeared = False
    motif_preserved = True

    def has_diaryl_ether(smiles):
        """Check if a molecule has a diaryl ether motif"""
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")
            return mol.HasSubstructMatch(diaryl_ether_pattern)
        return False

    def dfs_traverse(node, depth=0):
        nonlocal motif_appeared, motif_preserved

        if node["type"] == "mol" and node.get("smiles"):
            # Skip starting materials (in_stock)
            if node.get("in_stock", False):
                return

            has_motif = has_diaryl_ether(node["smiles"])

            # If motif has appeared before, check if it's preserved
            if motif_appeared and not has_motif:
                motif_preserved = False
                print(f"Diaryl ether motif lost in molecule: {node['smiles']}")

            # If this is the first time we see the motif, mark it
            if has_motif and not motif_appeared:
                motif_appeared = True
                print(f"Diaryl ether motif first appeared in: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the target molecule
    dfs_traverse(route)

    # If the motif never appeared, it can't be preserved
    if not motif_appeared:
        print("No diaryl ether motif found in any non-starting material molecule")
        return False

    return motif_preserved
