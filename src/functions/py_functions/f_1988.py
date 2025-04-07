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
    Detects if the synthesis route maintains a trifluoromethyl group throughout.
    """
    trifluoromethyl_present = False

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_present

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                trifluoromethyl_pattern = Chem.MolFromSmarts(
                    "[#6]-[#6]([#9])([#9])[#9]"
                )
                if mol.HasSubstructMatch(trifluoromethyl_pattern):
                    if depth == 0:  # Final product
                        trifluoromethyl_present = True
                        print("Trifluoromethyl group present in final product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return trifluoromethyl_present
