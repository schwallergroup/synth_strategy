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
    Detects the use of Weinreb amide intermediates in the synthesis.
    """
    weinreb_amide_found = False

    def dfs_traverse(node):
        nonlocal weinreb_amide_found

        if node["type"] == "mol" and "smiles" in node:
            # Check for Weinreb amide
            weinreb_pattern = Chem.MolFromSmarts("[#6][O][N]([#6])[C](=[O])[#6]")
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol and mol.HasSubstructMatch(weinreb_pattern):
                weinreb_amide_found = True
                print("Found Weinreb amide intermediate")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if weinreb_amide_found:
        print("Detected Weinreb amide intermediate strategy")
        return True
    return False
