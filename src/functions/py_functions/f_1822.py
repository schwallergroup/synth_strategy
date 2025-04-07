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
    Detects if the synthesis involves heterocyclic systems (thiazole, pyridine) throughout
    """
    thiazole_present = False
    pyridine_present = False

    def dfs_traverse(node, depth=0):
        nonlocal thiazole_present, pyridine_present

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for thiazole pattern
                thiazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#6][#16]1")
                if mol.HasSubstructMatch(thiazole_pattern):
                    thiazole_present = True
                    print(f"Found thiazole at depth {depth}")

                # Check for pyridine pattern
                pyridine_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#7]1")
                if mol.HasSubstructMatch(pyridine_pattern):
                    pyridine_present = True
                    print(f"Found pyridine at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both heterocycles are present
    return thiazole_present and pyridine_present
