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
    This function detects if the synthetic route follows a linear synthesis strategy
    without convergent steps.
    """
    # Track the maximum number of reactants in any step
    max_reactants = 0

    def dfs_traverse(node):
        nonlocal max_reactants

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count the number of reactants
            num_reactants = len(reactants)
            max_reactants = max(max_reactants, num_reactants)

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # If max_reactants is consistently low (≤2), it's likely a linear synthesis
    is_linear = max_reactants <= 2
    if is_linear:
        print(
            "Detected linear synthesis strategy with maximum",
            max_reactants,
            "reactants per step",
        )
    else:
        print(
            "Detected convergent synthesis with",
            max_reactants,
            "reactants in at least one step",
        )

    return is_linear
