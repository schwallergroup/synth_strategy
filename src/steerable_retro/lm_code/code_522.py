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
    Detects if the synthesis involves a heterocycle (like thiophene) that is maintained
    throughout the synthesis.
    """
    has_heterocycle = False

    def dfs_traverse(node):
        nonlocal has_heterocycle

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Check for thiophene pattern
                thiophene_pattern = Chem.MolFromSmarts("c1cscc1")
                if mol.HasSubstructMatch(thiophene_pattern):
                    has_heterocycle = True
                    print("Found thiophene heterocycle")

                # Check for other common heterocycles if needed
                # furan, pyrrole, etc.

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Heterocycle-containing synthesis detected: {has_heterocycle}")
    return has_heterocycle
