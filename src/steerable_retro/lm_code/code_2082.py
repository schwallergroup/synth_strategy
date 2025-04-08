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
    Detects if the synthesis involves preservation of a pre-existing stereocenter
    throughout the synthesis.
    """
    # Track if we found a stereocenter that's maintained
    stereocenter_preserved = False

    # Track molecules with stereocenters at each depth
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        nonlocal stereocenter_preserved

        if node["type"] == "mol":
            smiles = node["smiles"]
            if "@" in smiles:  # Simple check for stereochemistry in SMILES
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
                    if chiral_centers:
                        stereocenters_by_depth[depth] = chiral_centers
                        print(f"Found stereocenter at depth {depth}: {chiral_centers}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if stereocenters are preserved across multiple depths
    if len(stereocenters_by_depth) >= 2:
        # If we have stereocenters at multiple depths, consider them preserved
        stereocenter_preserved = True

    # Return True if stereocenters are preserved
    result = stereocenter_preserved
    print(f"Preservation of pre-existing stereocenter: {result}")
    return result
