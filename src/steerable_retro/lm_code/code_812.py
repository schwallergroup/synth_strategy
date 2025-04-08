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
    This function detects if the synthesis involves intermediates containing a nitrile group.
    """
    nitrile_found = False

    def dfs_traverse(node):
        nonlocal nitrile_found

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                # Nitrile SMARTS pattern
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                if mol and mol.HasSubstructMatch(nitrile_pattern):
                    nitrile_found = True
                    print("Nitrile group found in molecule")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return nitrile_found
