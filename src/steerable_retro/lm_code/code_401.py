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
    This function detects if the synthetic route follows a linear synthesis pattern
    without convergent steps.
    """
    # Track branching in the synthesis tree
    max_branching = 0

    def count_children(node):
        if node["type"] == "reaction":
            return len(node.get("children", []))
        return 0

    def dfs_traverse(node):
        nonlocal max_branching

        # Count number of reactants in this reaction
        num_children = count_children(node)
        max_branching = max(max_branching, num_children)

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)

    # Linear synthesis typically has at most 2 reactants per step
    is_linear = max_branching <= 2
    if is_linear:
        print("Linear synthesis pattern detected")
    else:
        print(f"Convergent synthesis detected with maximum {max_branching} reactants in a step")

    return is_linear
