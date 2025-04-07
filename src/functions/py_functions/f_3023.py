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
    Detects if the synthesis route preserves stereochemistry throughout the synthesis.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Count chiral centers with defined stereochemistry
                chiral_centers = 0
                for atom in mol.GetAtoms():
                    if atom.HasProp("_CIPCode"):  # R or S configuration
                        chiral_centers += 1

                if chiral_centers > 0:
                    stereocenters_by_depth[depth] = chiral_centers
                    print(f"Found {chiral_centers} stereocenters at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if stereocenters are preserved throughout (present at multiple depths)
    preserved = len(stereocenters_by_depth) >= 2

    if preserved:
        print(
            f"Stereocenters preserved across multiple steps: {stereocenters_by_depth}"
        )

    return preserved
