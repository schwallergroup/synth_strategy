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
    This function detects if the synthesis route involves a tetralin scaffold
    that is maintained throughout the synthesis.
    """
    tetralin_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal tetralin_count

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Tetralin scaffold pattern
                tetralin_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6][#6][#6]2[#6][#6][#6][#6][#6]2[#6]1"
                )
                if mol.HasSubstructMatch(tetralin_pattern):
                    tetralin_count += 1
                    print(f"Detected tetralin scaffold at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # Return True if tetralin scaffold appears in multiple intermediates
    return tetralin_count >= 2
