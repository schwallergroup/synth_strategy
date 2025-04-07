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
    This function detects if a piperazine scaffold is maintained throughout
    the synthesis.
    """
    # Track piperazine presence at different depths
    piperazine_depths = []

    # SMARTS pattern for piperazine
    piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(piperazine_pattern):
                piperazine_depths.append(depth)
                print(f"Piperazine detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if piperazine is maintained through multiple depths
    if len(set(piperazine_depths)) >= 3:  # Present in at least 3 different depths
        print(f"Piperazine scaffold maintained through multiple steps: {sorted(piperazine_depths)}")
        return True
    return False
