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
    This function detects a synthetic strategy involving a thioether-containing linker.
    """
    found_thioether_linker = False

    def dfs_traverse(node):
        nonlocal found_thioether_linker

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None:
                # Look for thioether pattern with carbon chains on both sides
                thioether_pattern = Chem.MolFromSmarts("[#6][S][#6][#6]")
                if len(mol.GetSubstructMatches(thioether_pattern)) > 0:
                    found_thioether_linker = True
                    print(f"Detected thioether-containing linker in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Thioether-containing linker detected: {found_thioether_linker}")
    return found_thioether_linker
