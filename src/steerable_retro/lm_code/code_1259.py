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
    Detects if the synthesis follows a linear strategy without convergent steps.
    """
    # Track the maximum branching factor
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count number of reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            num_reactants = len([r for r in reactants if r])

            # Update max branching
            max_branching = max(max_branching, num_reactants)
            print(f"Reaction has {num_reactants} reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis has max branching factor of 1
    return max_branching <= 2  # Allow for solvents/reagents
