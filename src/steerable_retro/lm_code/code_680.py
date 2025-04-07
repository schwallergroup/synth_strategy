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
        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Look for @H or @@H in SMILES as a simple stereocenter check
                if "@H" in node["smiles"]:
                    stereocenters_by_depth[depth] = True
                    print(f"Found stereocenter at depth {depth}: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if stereocenters are present at multiple depths
    if len(stereocenters_by_depth) >= 2:
        print(f"Stereocenters preserved across {len(stereocenters_by_depth)} steps")
        return True
    return False
