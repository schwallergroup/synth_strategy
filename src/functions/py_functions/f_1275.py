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
    This function detects if stereochemistry is preserved throughout the synthesis
    without creation of new stereocenters.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        nonlocal stereocenters_by_depth

        if node["type"] == "mol":
            # Count stereocenters in the molecule
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                stereocenters_by_depth[depth] = len(chiral_centers)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if stereocenters are preserved (no new ones created)
    if not stereocenters_by_depth:
        return False

    # Get the number of stereocenters at the final product (depth 0)
    final_stereocenters = stereocenters_by_depth.get(0, 0)

    # Check if all intermediates have the same or fewer stereocenters
    for depth, count in stereocenters_by_depth.items():
        if depth > 0 and count > final_stereocenters:
            print(f"New stereocenters created at depth {depth}")
            return False

    print("Stereochemistry preserved throughout synthesis")
    return True
