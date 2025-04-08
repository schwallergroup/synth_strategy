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
    This function detects if a synthesis preserves an amide functional group
    throughout the synthesis.
    """
    # SMARTS pattern for amide group
    amide_pattern = Chem.MolFromSmarts("[NH][C](=O)")

    # Track molecules at each depth
    molecules_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            # Store molecule at this depth
            if depth not in molecules_by_depth:
                molecules_by_depth[depth] = []
            molecules_by_depth[depth].append(node["smiles"])

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if amide group is present at all depths
    amide_preserved = True
    for depth, molecules in molecules_by_depth.items():
        depth_has_amide = False
        for smiles in molecules:
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol and mol.HasSubstructMatch(amide_pattern):
                    depth_has_amide = True
                    break
            except:
                continue

        if (
            not depth_has_amide and molecules
        ):  # If we have molecules at this depth but none with amide
            amide_preserved = False
            break

    if amide_preserved:
        print("Amide functional group is preserved throughout synthesis")
    else:
        print("Amide functional group is not preserved throughout synthesis")

    return amide_preserved
