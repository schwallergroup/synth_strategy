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
    Detects if the synthesis follows a linear strategy with sequential functionalization
    rather than a convergent approach.
    """
    reaction_count = 0
    max_reactants_per_step = 0

    def dfs_traverse(node):
        nonlocal reaction_count, max_reactants_per_step

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count number of reactants
            num_reactants = len(reactants)
            max_reactants_per_step = max(max_reactants_per_step, num_reactants)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Linear synthesis typically has few reactants per step (mostly 1-2)
    # and follows a sequential pattern
    is_linear = reaction_count >= 3 and max_reactants_per_step <= 2

    if is_linear:
        print(f"Detected linear synthesis strategy with {reaction_count} steps")

    return is_linear
