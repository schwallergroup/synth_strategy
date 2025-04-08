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
    This function detects a strategy where heterocycles (tetrazole and piperidine)
    are preserved throughout the synthesis.
    """
    tetrazole_depths = []
    piperidine_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for tetrazole
                    tetrazole_pattern = Chem.MolFromSmarts("[n]1[n][n][n][c]1")
                    if mol.HasSubstructMatch(tetrazole_pattern):
                        tetrazole_depths.append(depth)
                        print(f"Tetrazole found at depth {depth}")

                    # Check for piperidine
                    piperidine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#6][#6][#6]1")
                    if mol.HasSubstructMatch(piperidine_pattern):
                        piperidine_depths.append(depth)
                        print(f"Piperidine found at depth {depth}")
            except:
                print(f"Error processing molecule SMILES at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if both heterocycles are present in the final product (depth 0)
    # and preserved from earlier steps
    if (
        0 in tetrazole_depths
        and len(tetrazole_depths) > 1
        and 0 in piperidine_depths
        and len(piperidine_depths) > 1
    ):
        print("Heterocycle preservation strategy detected")
        return True
    return False
