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
    This function detects a linear synthesis strategy where each step
    has predominantly two reactants (one main substrate and one reagent).
    """
    linear_steps = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal linear_steps, total_steps

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            total_steps += 1
            # A linear step typically has 1-2 reactants
            if len(reactants) <= 2:
                linear_steps += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present (at least 80% of steps are linear)
    strategy_present = total_steps > 0 and (linear_steps / total_steps) >= 0.8

    print(
        f"Linear synthesis strategy detected: {strategy_present} ({linear_steps}/{total_steps} linear steps)"
    )
    return strategy_present
