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
    This function detects if the synthetic route involves compounds containing a piperazine ring.
    """
    piperazine_found = False

    def dfs_traverse(node):
        nonlocal piperazine_found

        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for piperazine ring
                    piperazine_pattern = Chem.MolFromSmarts("[N]1[CH2][CH2][N][CH2][CH2]1")
                    if mol.HasSubstructMatch(piperazine_pattern):
                        print("Piperazine-containing compound detected")
                        piperazine_found = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return piperazine_found
