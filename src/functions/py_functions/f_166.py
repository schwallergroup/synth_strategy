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
    This function detects a synthetic strategy involving linear assembly of fragments
    rather than convergent synthesis.
    """
    # Track reaction depths and number of fragments at each step
    reaction_data = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Record number of fragments at this depth
                reaction_data.append((depth, len(reactants_smiles)))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze reaction data to determine if synthesis is linear
    # In linear synthesis, most steps have 2 or fewer reactants
    if not reaction_data:
        return False

    linear_steps = sum(1 for _, num_reactants in reaction_data if num_reactants <= 2)
    total_steps = len(reaction_data)

    # If more than 75% of steps are linear (2 or fewer reactants), consider it a linear strategy
    is_linear = (linear_steps / total_steps) >= 0.75

    if is_linear:
        print(
            f"Linear fragment assembly detected: {linear_steps}/{total_steps} steps are linear"
        )

    return is_linear
