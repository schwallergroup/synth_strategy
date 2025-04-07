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
    This function detects preservation of diaryl ether linkage
    throughout the synthetic route.
    """
    # Track if we find the diaryl ether in the final product
    diaryl_ether_preserved = False

    def dfs_traverse(node):
        nonlocal diaryl_ether_preserved

        if node["type"] == "mol" and node.get("in_stock", False) == False:
            # This is likely the final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for diaryl ether
                diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")
                if mol.HasSubstructMatch(diaryl_ether_pattern):
                    print("Found diaryl ether in final product")
                    diaryl_ether_preserved = True

        # Check if diaryl ether is present in starting materials
        if node["type"] == "mol" and node.get("in_stock", False) == True:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for diaryl ether in starting materials
                diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")
                if mol.HasSubstructMatch(diaryl_ether_pattern):
                    print("Found diaryl ether in starting material")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if diaryl ether is preserved
    return diaryl_ether_preserved
