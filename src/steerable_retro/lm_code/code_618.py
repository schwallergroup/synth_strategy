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
    This function detects a linear build-up strategy where fragments are added sequentially
    rather than in a convergent manner.
    """
    reaction_count = 0
    max_reactants_per_step = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, max_reactants_per_step

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Count number of reactants
            num_reactants = len(reactants)
            max_reactants_per_step = max(max_reactants_per_step, num_reactants)

            print(f"Reaction at depth {depth} has {num_reactants} reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If most reactions have 2 or fewer reactants and we have multiple steps,
    # this suggests a linear build-up rather than convergent synthesis
    if reaction_count >= 3 and max_reactants_per_step <= 2:
        print(f"Detected linear build-up strategy with {reaction_count} reactions")
        return True

    return False
