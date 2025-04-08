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
    This function detects if the synthesis involves a pyrazole core with a
    trifluoromethoxy-containing aromatic fragment throughout the synthesis.
    """
    # Define SMARTS patterns
    pyrazole_pattern = "[n]1[n][c][c][c]1"
    trifluoromethoxy_pattern = "[O][C]([F])([F])[F]"

    pyrazole_present = False
    trifluoromethoxy_present = False

    def dfs_traverse(node, depth=0):
        nonlocal pyrazole_present, trifluoromethoxy_present

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for pyrazole
                    pyrazole_matcher = Chem.MolFromSmarts(pyrazole_pattern)
                    if mol.HasSubstructMatch(pyrazole_matcher):
                        pyrazole_present = True
                        print(f"Found pyrazole at depth {depth}")

                    # Check for trifluoromethoxy
                    trifluoromethoxy_matcher = Chem.MolFromSmarts(trifluoromethoxy_pattern)
                    if mol.HasSubstructMatch(trifluoromethoxy_matcher):
                        trifluoromethoxy_present = True
                        print(f"Found trifluoromethoxy at depth {depth}")
            except Exception as e:
                print(f"Error analyzing molecule at depth {depth}: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True if both substructures are present
    result = pyrazole_present and trifluoromethoxy_present
    if result:
        print("Detected pyrazole and trifluoromethoxy-containing synthesis")
    return result
