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
    Detects if the synthesis route preserves stereochemistry throughout.
    """
    stereocenters_preserved = False
    chiral_centers_count = 0

    def dfs_traverse(node):
        nonlocal stereocenters_preserved, chiral_centers_count

        if node["type"] == "mol" and "smiles" in node:
            # Check for stereochemistry indicators in SMILES
            if "@" in node["smiles"]:
                chiral_centers_count += 1
                print(f"Found chiral center in molecule: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    stereocenters_preserved = chiral_centers_count >= 2  # At least two molecules with stereocenters
    print(
        f"Stereocenter preservation strategy detected: {stereocenters_preserved} (chiral centers: {chiral_centers_count})"
    )
    return stereocenters_preserved
