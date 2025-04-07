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
    Detects if the synthesis follows a linear strategy (each step builds upon the previous)
    rather than a convergent strategy.
    """
    # Track the maximum branching factor in the synthetic tree
    max_branching = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count the number of non-trivial reactants (molecules that are synthesized, not starting materials)
            non_trivial_reactants = 0

            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_trivial_reactants += 1

            max_branching = max(max_branching, non_trivial_reactants)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # If max_branching is 1, the synthesis is linear
    is_linear = max_branching <= 1

    if is_linear:
        print("Detected linear synthesis strategy (no convergent steps)")
    else:
        print(f"Detected convergent synthesis with maximum branching factor of {max_branching}")

    return is_linear
