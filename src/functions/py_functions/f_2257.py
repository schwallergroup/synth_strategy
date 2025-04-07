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
    This function detects the use of trifluoromethyl-containing building blocks
    in the synthesis route.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if node["type"] == "mol" and node.get("in_stock", False):
            # Check if this starting material contains a trifluoromethyl group
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol is not None:
                    trifluoromethyl_pattern = Chem.MolFromSmarts("[#6][C]([F])([F])[F]")
                    if mol.HasSubstructMatch(trifluoromethyl_pattern):
                        print(
                            f"Found trifluoromethyl-containing building block at depth {depth}"
                        )
                        found_pattern = True
            except Exception as e:
                print(f"Error processing molecule SMILES: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_pattern
