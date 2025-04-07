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
    This function detects if aromatic halides (Cl, Br) are preserved
    throughout the synthesis.
    """
    has_final_halides = False
    has_initial_halides = False

    def dfs_traverse(node, depth=0, is_leaf=False):
        nonlocal has_final_halides, has_initial_halides

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for aromatic halides
                ar_cl_pattern = Chem.MolFromSmarts("c-[#17]")
                ar_br_pattern = Chem.MolFromSmarts("c-[#35]")

                has_halides = mol.HasSubstructMatch(ar_cl_pattern) or mol.HasSubstructMatch(
                    ar_br_pattern
                )

                if depth == 0:  # Final product
                    has_final_halides = has_halides
                    print(f"Final product has halides: {has_halides}")

                if node.get("in_stock", False) or is_leaf:  # Starting material
                    if has_halides:
                        has_initial_halides = True
                        print("Starting material has halides")

        # Traverse children
        children = node.get("children", [])
        is_leaf = len(children) == 0
        for child in children:
            dfs_traverse(child, depth + 1, is_leaf)

    # Start traversal from root
    dfs_traverse(route)
    return has_final_halides and has_initial_halides
