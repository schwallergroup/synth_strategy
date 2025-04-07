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
    Detects synthesis involving methoxy-substituted aromatic compounds.
    """
    # Track if we found the pattern
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        # Check for methoxy aromatic pattern in molecules
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    methoxy_pattern = Chem.MolFromSmarts("[c]-[O]-[CH3]")
                    if mol.HasSubstructMatch(methoxy_pattern):
                        found_pattern = True
                        print(f"Found methoxy-substituted aromatic at depth {depth}")
            except:
                pass

        # Traverse children
        if "children" in node:
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Methoxy aromatic containing synthesis: {found_pattern}")
    return found_pattern
