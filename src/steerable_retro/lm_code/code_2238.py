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
    Detects if the synthesis involves a tetrazole ring that is maintained throughout
    """
    has_tetrazole_in_final = False
    has_tetrazole_in_intermediate = False

    def dfs_traverse(node, depth=0):
        nonlocal has_tetrazole_in_final, has_tetrazole_in_intermediate

        if node["type"] == "mol":
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Tetrazole pattern
                    tetrazole_pattern = Chem.MolFromSmarts("c1nnn[nH]1")
                    tetrazole_pattern2 = Chem.MolFromSmarts("c1nnn[n]1-[#6]")

                    if mol.HasSubstructMatch(tetrazole_pattern) or mol.HasSubstructMatch(
                        tetrazole_pattern2
                    ):
                        if depth == 0:  # Final product
                            has_tetrazole_in_final = True
                            print("Final product contains tetrazole")
                        else:  # Intermediate
                            has_tetrazole_in_intermediate = True
                            print(f"Intermediate at depth {depth} contains tetrazole")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_tetrazole_in_final and has_tetrazole_in_intermediate
