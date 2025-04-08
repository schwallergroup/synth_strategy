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
    This function detects if the final product contains an indole scaffold.
    """
    contains_indole = False

    def dfs_traverse(node):
        nonlocal contains_indole

        if node["type"] == "mol" and "smiles" in node and not node.get("children"):
            # This is a leaf node (final product)
            indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(indole_pattern):
                contains_indole = True
                print("Final product contains indole scaffold")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return contains_indole
