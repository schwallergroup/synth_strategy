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
    Detects synthesis routes involving isoxazole heterocycles.
    """
    isoxazole_detected = False

    def dfs_traverse(node):
        nonlocal isoxazole_detected

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is not None:
                        # SMARTS pattern for isoxazole
                        isoxazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#8][#6][#6]1")
                        if mol.HasSubstructMatch(isoxazole_pattern):
                            isoxazole_detected = True
                            print(f"Detected isoxazole in molecule: {smiles}")
                except Exception as e:
                    print(f"Error in isoxazole detection: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return isoxazole_detected
