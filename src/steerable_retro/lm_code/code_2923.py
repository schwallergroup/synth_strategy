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
    This function detects if halogens (F, I) present in starting materials
    are retained in the final product.
    """
    starting_halogens = set()
    final_halogens = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for halogens
                    f_pattern = Chem.MolFromSmarts("[F]")
                    i_pattern = Chem.MolFromSmarts("[I]")

                    has_f = mol.HasSubstructMatch(f_pattern)
                    has_i = mol.HasSubstructMatch(i_pattern)

                    if node.get("in_stock", False):  # Starting material
                        if has_f:
                            starting_halogens.add("F")
                        if has_i:
                            starting_halogens.add("I")
                    elif depth == 0:  # Final product
                        if has_f:
                            final_halogens.add("F")
                        if has_i:
                            final_halogens.add("I")
            except Exception as e:
                print(f"Error in halogen analysis: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if all starting halogens are retained in final product
    if starting_halogens and starting_halogens.issubset(final_halogens):
        print(f"Detected halogen retention: {starting_halogens} retained in final product")
        return True
    return False
