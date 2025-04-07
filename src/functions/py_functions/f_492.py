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
    Detects if a stereocenter is preserved throughout the synthesis.
    """
    stereocenters_by_depth = {}

    def dfs_traverse(node, depth=0):
        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]

            # Check for stereochemistry markers in SMILES
            if "@" in product_smiles:
                stereocenters_by_depth[depth] = True
                print(f"Stereocenter found at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if stereocenters are present in at least 2 consecutive depths
    preserved = len(stereocenters_by_depth) >= 2

    if preserved:
        print(f"Stereocenter preserved through {len(stereocenters_by_depth)} reactions")

    return preserved
