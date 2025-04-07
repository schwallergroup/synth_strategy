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
    This function detects if the synthesis route uses azide as an intermediate.
    """
    azide_found = False

    def dfs_traverse(node):
        nonlocal azide_found

        if node["type"] == "mol":
            # Check if molecule contains azide group
            if node["smiles"] and "[N-]=[N+]=[N" in node["smiles"]:
                azide_found = True
                print(f"Azide intermediate found: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    return azide_found
