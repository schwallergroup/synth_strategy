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
    Detects if the synthesis route involves compounds containing trifluoromethyl groups.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "mol":
            smiles = node["smiles"]
            if "C(F)(F)F" in smiles or "CF3" in smiles:
                print("Detected trifluoromethyl group")
                result = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return result
