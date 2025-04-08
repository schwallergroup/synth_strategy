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
    This function detects a linear synthesis strategy as opposed to convergent.
    It analyzes the structure of the synthetic tree to determine if it's primarily linear.
    """
    max_branching = 0

    def count_reactants(rsmi):
        reactants = rsmi.split(">")[0].split(".")
        return len([r for r in reactants if r])

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                num_reactants = count_reactants(rsmi)
                max_branching = max(max_branching, num_reactants)
                print(f"Reaction with {num_reactants} reactants found")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If max branching is <= 2, consider it a linear synthesis
    # (allowing for one main substrate and one reagent)
    return max_branching <= 2
