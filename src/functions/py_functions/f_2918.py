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
    Detects if the synthetic route includes a trifluoromethyl group.
    """
    trifluoromethyl_found = False

    def dfs_traverse(node):
        nonlocal trifluoromethyl_found

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[#6]([F])([F])[F]")
                ):
                    print("Trifluoromethyl group detected")
                    trifluoromethyl_found = True
            except:
                print("Error processing molecule SMILES for trifluoromethyl detection")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return trifluoromethyl_found
