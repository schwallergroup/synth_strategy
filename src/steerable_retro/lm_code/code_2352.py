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
    This function detects a synthetic strategy involving azide intermediates.
    """
    azide_intermediate_found = False

    def dfs_traverse(node):
        nonlocal azide_intermediate_found

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for azide group
                azide_pattern = Chem.MolFromSmarts("[N-]=[N+]=N")
                if mol.HasSubstructMatch(azide_pattern):
                    print("Found azide intermediate")
                    azide_intermediate_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return azide_intermediate_found
