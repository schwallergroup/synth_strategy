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
    Detects a synthetic strategy involving multiple halogen-containing intermediates
    with different halogens (Br, Cl, I) used as synthetic handles.
    """
    halogen_types = set()

    def dfs_traverse(node, depth=0):
        nonlocal halogen_types

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]

            # Check for different halogens
            if "Br" in smiles:
                halogen_types.add("Br")
            if "Cl" in smiles:
                halogen_types.add("Cl")
            if "I" in smiles:
                halogen_types.add("I")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Found {len(halogen_types)} different halogen types: {halogen_types}")
    return len(halogen_types) >= 2
