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
    This function detects if a Stille coupling reaction is used in the synthesis.
    """
    found_stille = False

    def dfs_traverse(node, depth=0):
        nonlocal found_stille

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Stille coupling (aryl-I + vinyl-Sn → aryl-vinyl)
            if any("[Sn]" in r for r in reactants) and any("I" in r for r in reactants):
                print(f"Found Stille coupling at depth {depth}")
                found_stille = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if found_stille:
        print("Stille coupling strategy detected")
    else:
        print("No Stille coupling detected")

    return found_stille
