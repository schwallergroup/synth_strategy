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
    This function detects the presence of methoxy-substituted and
    fluorinated aromatic compounds in the synthesis.
    """
    # SMARTS patterns
    methoxy_aromatic = Chem.MolFromSmarts("[OX2]([CH3])[c]")
    fluoro_aromatic = Chem.MolFromSmarts("[F][c]")

    # Track if we found each pattern
    found_methoxy = False
    found_fluoro = False

    def dfs_traverse(node):
        nonlocal found_methoxy, found_fluoro

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if mol.HasSubstructMatch(methoxy_aromatic):
                    found_methoxy = True
                    print("Found methoxy-substituted aromatic")
                if mol.HasSubstructMatch(fluoro_aromatic):
                    found_fluoro = True
                    print("Found fluorinated aromatic")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found both patterns
    result = found_methoxy and found_fluoro
    print(f"Aromatic methoxy and fluorine pattern: {result}")
    return result
