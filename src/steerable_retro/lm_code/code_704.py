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
    Detects a synthetic strategy involving a difluorinated aromatic system
    that is maintained throughout the synthesis.
    """
    difluoro_aromatic_present = False

    def dfs_traverse(node, depth=0):
        nonlocal difluoro_aromatic_present

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if not smiles:
                return

            # Check for difluorinated aromatic system
            # Pattern looks for an aromatic ring with two fluorine atoms
            difluoro_pattern = Chem.MolFromSmarts("c([F])cc([F])")

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(difluoro_pattern):
                    difluoro_aromatic_present = True
                    print(f"Found difluorinated aromatic system at depth {depth}")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return difluoro_aromatic_present
