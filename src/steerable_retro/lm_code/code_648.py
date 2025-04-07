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
    This function detects preservation of stereochemistry throughout the synthesis.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Find chiral centers
                chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
                if chiral_centers:
                    print(f"Found {len(chiral_centers)} chiral centers at depth {depth}")
                    stereocenters_by_depth[depth] = len(chiral_centers)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if stereocenter count is consistent throughout
    if len(stereocenters_by_depth) >= 2:  # At least two steps with stereocenters
        stereocenter_counts = list(stereocenters_by_depth.values())
        if all(count == stereocenter_counts[0] for count in stereocenter_counts):
            print("Stereochemistry preserved throughout synthesis")
            return True

    return False
