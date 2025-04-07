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
    Detects a strategy where a morpholine amide functional group is retained
    throughout multiple steps of the synthesis.
    """
    morpholine_amide_depths = []

    # SMARTS pattern for morpholine amide
    morpholine_amide_pattern = Chem.MolFromSmarts(
        "[#6](=[O])[#7]1[#6][#6][#8][#6][#6]1"
    )

    def dfs_traverse(node):
        if node["type"] == "mol":
            if node["smiles"]:
                mol = Chem.MolFromSmiles(node["smiles"])
                depth = node.get("depth", 0)

                if mol and mol.HasSubstructMatch(morpholine_amide_pattern):
                    morpholine_amide_depths.append(depth)
                    print(f"Found morpholine amide at depth {depth}: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if morpholine amide is present at multiple depths
    if len(morpholine_amide_depths) >= 2:
        # Sort depths to check range
        morpholine_amide_depths.sort()

        # Check if the morpholine amide spans a significant portion of the synthesis
        depth_range = morpholine_amide_depths[-1] - morpholine_amide_depths[0]
        if depth_range >= 2:  # Present across at least 3 steps
            print(f"Morpholine amide retained across {depth_range+1} steps")
            return True

    return False
