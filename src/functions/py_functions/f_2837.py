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
    This function detects a linear synthesis strategy (no convergent steps).
    """
    # Track the maximum branching factor
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count number of mol children (reactants)
            mol_children = sum(
                1 for child in node.get("children", []) if child["type"] == "mol"
            )
            max_branching = max(max_branching, mol_children)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Linear synthesis has max branching factor of 2 (binary reactions)
    is_linear = max_branching <= 2
    if is_linear:
        print("Found linear synthesis strategy")

    return is_linear
