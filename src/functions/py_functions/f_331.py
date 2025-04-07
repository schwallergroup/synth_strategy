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
    This function detects if the synthesis involves a morpholine-containing scaffold.
    """
    morpholine_found = False

    def dfs_traverse(node):
        nonlocal morpholine_found

        if node["type"] == "mol":
            # Check if molecule contains morpholine
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(
                Chem.MolFromSmarts("[N]1[C][C][O][C][C]1")
            ):
                print("Morpholine-containing scaffold detected")
                morpholine_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return morpholine_found
