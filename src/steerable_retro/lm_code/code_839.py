#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects if the synthesis involves heterocyclic compounds,
    particularly pyrazole and pyridine rings.
    """
    has_pyrazole = False
    has_pyridine = False

    def dfs_traverse(node):
        nonlocal has_pyrazole, has_pyridine

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Check for pyrazole
                    pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c,n][c,n][c,n]1")
                    if mol.HasSubstructMatch(pyrazole_pattern):
                        print("Found pyrazole heterocycle")
                        has_pyrazole = True

                    # Check for pyridine
                    pyridine_pattern = Chem.MolFromSmarts("[n]1[c][c][c][c][c]1")
                    if mol.HasSubstructMatch(pyridine_pattern):
                        print("Found pyridine heterocycle")
                        has_pyridine = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_pyrazole or has_pyridine  # Either heterocycle is present
