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
    Detects if the synthesis route involves multiple heterocycles
    (indole, pyrazolopyrimidine, morpholine).
    """
    heterocycle_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_count

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if not mol:
                return

            # Define patterns for different heterocycles
            indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
            pyrazolopyrimidine_pattern = Chem.MolFromSmarts("c1nc2cncnc2[nH]1")
            morpholine_pattern = Chem.MolFromSmarts("[#8]1[#6][#6][#7][#6][#6]1")

            # Check for each heterocycle
            if mol.HasSubstructMatch(indole_pattern):
                heterocycle_count += 1
                print("Detected indole heterocycle")

            if mol.HasSubstructMatch(pyrazolopyrimidine_pattern):
                heterocycle_count += 1
                print("Detected pyrazolopyrimidine heterocycle")

            if mol.HasSubstructMatch(morpholine_pattern):
                heterocycle_count += 1
                print("Detected morpholine heterocycle")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return heterocycle_count >= 2  # Return True if at least 2 different heterocycles are present
