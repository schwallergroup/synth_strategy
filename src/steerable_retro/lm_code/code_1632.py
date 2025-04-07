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
    Detects if the synthesis route involves multiple convergent steps.
    """
    convergent_steps = 0

    def dfs_traverse(node):
        nonlocal convergent_steps

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If there are multiple reactants, it's a convergent step
            if len(reactants) > 1:
                convergent_steps += 1
                print(f"Convergent step detected, total: {convergent_steps}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy requires at least 2 convergent steps
    result = convergent_steps >= 2
    print(f"Multiple convergent steps strategy detected: {result}")
    return result
