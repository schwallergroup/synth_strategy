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
    This function detects if the synthesis uses a thiophene-containing building block.
    """
    contains_thiophene = False

    def dfs_traverse(node):
        nonlocal contains_thiophene

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    thiophene_pattern = Chem.MolFromSmarts("c1cscc1")
                    if mol.HasSubstructMatch(thiophene_pattern):
                        contains_thiophene = True
                        print("Detected thiophene-containing building block")
            except:
                print("Error processing molecule SMILES in thiophene detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return contains_thiophene
