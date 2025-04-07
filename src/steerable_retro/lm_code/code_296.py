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
    Detects synthesis routes that incorporate a morpholine group.
    """
    contains_morpholine = False

    def dfs_traverse(node):
        nonlocal contains_morpholine

        if node["type"] == "mol" and node.get("smiles"):
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Morpholine pattern
                    morpholine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#8][#6][#6]1")
                    if mol.HasSubstructMatch(morpholine_pattern):
                        print("Found morpholine-containing molecule")
                        contains_morpholine = True
            except:
                print("Error processing molecule SMILES for morpholine detection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return contains_morpholine
