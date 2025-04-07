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
    This function detects the use of nitrile as a key intermediate in the synthesis.
    """
    nitrile_present = False

    def dfs_traverse(node):
        nonlocal nitrile_present

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            # Nitrile pattern
            nitrile_pattern = Chem.MolFromSmarts("[C;X2]#[N;X1]")

            if mol and mol.HasSubstructMatch(nitrile_pattern):
                nitrile_present = True
                print(f"Nitrile intermediate detected: {smiles}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitrile_present
