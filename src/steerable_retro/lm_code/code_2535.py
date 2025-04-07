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
    This function detects if a nitrile functional group is maintained through
    multiple steps of the synthesis.
    """
    # Track nitrile presence at different depths
    nitrile_depths = []

    # SMARTS pattern for nitrile group
    nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(nitrile_pattern):
                nitrile_depths.append(depth)
                print(f"Nitrile detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if nitrile is maintained through multiple depths
    if len(set(nitrile_depths)) >= 3:  # Present in at least 3 different depths
        print(f"Nitrile maintained through multiple steps: {sorted(nitrile_depths)}")
        return True
    return False
