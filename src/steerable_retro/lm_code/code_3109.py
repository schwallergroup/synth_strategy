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
    Detects if the synthesis involves a trifluoromethyl group.
    """
    has_trifluoromethyl = False

    def dfs_traverse(node):
        nonlocal has_trifluoromethyl

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if not mol:
                return

            # Check for trifluoromethyl pattern
            cf3_patt = Chem.MolFromSmarts("[C]([F])([F])[F]")
            if mol.HasSubstructMatch(cf3_patt):
                has_trifluoromethyl = True
                print("Detected trifluoromethyl group")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return has_trifluoromethyl
