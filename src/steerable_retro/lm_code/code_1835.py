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
    Detects if the synthesis involves a nitrile-containing aromatic fragment
    that is maintained throughout the synthesis.
    """
    # Track if we found a nitrile-containing fragment
    found_nitrile = False
    # Track at what depths we find the nitrile
    nitrile_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal found_nitrile, nitrile_depths

        if node["type"] == "mol":
            if "smiles" in node:
                # Check for aromatic ring with nitrile
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    nitrile_pattern = Chem.MolFromSmarts("c:c-[#6]#[#7]")
                    if mol.HasSubstructMatch(nitrile_pattern):
                        found_nitrile = True
                        nitrile_depths.append(depth)
                        print(f"Found nitrile-containing aromatic fragment at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if nitrile is present throughout the synthesis
    # We should find it at multiple depths including depth 0 (final product)
    if found_nitrile and 0 in nitrile_depths and len(nitrile_depths) > 1:
        print("Found nitrile-containing aromatic fragment maintained throughout synthesis")
        return True

    return False
