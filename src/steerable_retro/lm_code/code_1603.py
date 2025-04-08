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
    This function detects a strategy using halogen-containing building blocks
    (specifically fluorophenyl) that are carried through the synthesis.
    """
    fluorophenyl_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for fluorophenyl group
                    fluorophenyl_pattern = Chem.MolFromSmarts("[c]1[c][c][c]([F])[c][c]1")
                    if mol.HasSubstructMatch(fluorophenyl_pattern):
                        fluorophenyl_depths.append(depth)
                        print(f"Fluorophenyl group found at depth {depth}")
            except:
                print(f"Error processing molecule SMILES at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if fluorophenyl is present in final product and in early steps
    if 0 in fluorophenyl_depths and max(fluorophenyl_depths) >= 3:
        print("Halogen building block strategy detected")
        return True
    return False
