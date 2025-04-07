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
    Detects if the synthesis route involves multiple nitrogen-containing
    heterocycles (pyrimidine, pyrazole, etc.).
    """
    heterocycle_count = 0

    def dfs_traverse(node):
        nonlocal heterocycle_count

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Check for pyrimidine
                    pyrimidine_pattern = Chem.MolFromSmarts("n1cnccc1")
                    # Check for pyrazole
                    pyrazole_pattern = Chem.MolFromSmarts("n1ncc[c,n]1")
                    # Check for piperidine
                    piperidine_pattern = Chem.MolFromSmarts("N1CCCCC1")

                    unique_heterocycles = 0
                    if mol.HasSubstructMatch(pyrimidine_pattern):
                        unique_heterocycles += 1
                    if mol.HasSubstructMatch(pyrazole_pattern):
                        unique_heterocycles += 1
                    if mol.HasSubstructMatch(piperidine_pattern):
                        unique_heterocycles += 1

                    heterocycle_count = max(heterocycle_count, unique_heterocycles)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_count >= 2  # Return True if at least 2 different heterocycles detected
