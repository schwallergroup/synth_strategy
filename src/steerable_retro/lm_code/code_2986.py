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
    This function detects if the synthesis involves an indole-containing compound.
    """
    has_indole = False

    def dfs_traverse(node):
        nonlocal has_indole

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Indole pattern
                indole_pattern = Chem.MolFromSmarts("c1cccc2c1nc[c]2")
                if mol.HasSubstructMatch(indole_pattern):
                    has_indole = True
                    print("Indole scaffold detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_indole
