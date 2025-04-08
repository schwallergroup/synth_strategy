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
    This function detects if the synthesis follows a linear strategy without convergent steps.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")

            # If there are more than 2 reactants, it's likely a convergent step
            if len(reactants) > 2:
                is_linear = False
                print(f"Found convergent step with {len(reactants)} reactants at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Linear synthesis strategy: {is_linear}")
    return is_linear
