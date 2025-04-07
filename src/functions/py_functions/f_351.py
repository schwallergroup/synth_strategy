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
    Detects if the synthesis route involves sequential elaboration of heterocyclic systems,
    specifically thiophene and pyrazole rings.
    """
    # Track presence of heterocycles
    thiophene_present = False
    pyrazole_present = False

    def dfs_traverse(node, depth=0):
        nonlocal thiophene_present, pyrazole_present

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if not smiles:
                return

            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return

            # Check for thiophene
            thiophene_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#16][#6]1")
            if mol.HasSubstructMatch(thiophene_pattern):
                thiophene_present = True
                print("Detected thiophene ring")

            # Check for pyrazole
            pyrazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#7][#6][#6]1")
            if mol.HasSubstructMatch(pyrazole_pattern):
                pyrazole_present = True
                print("Detected pyrazole ring")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return true if both heterocycles are present
    return thiophene_present and pyrazole_present
