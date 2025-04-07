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
    This function detects a linear synthesis strategy with sequential transformations
    and no convergent steps.
    """
    # Track the maximum number of reactants in any step
    max_reactants = 1

    def dfs_traverse(node):
        nonlocal max_reactants

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # Count number of reactants
            num_reactants = len(reactants)
            max_reactants = max(max_reactants, num_reactants)

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If max_reactants is consistently low (≤2), it's likely a linear synthesis
    # Most reactions have at least 2 reactants (main substrate + reagent)
    print(f"Maximum reactants in any step: {max_reactants}")
    return max_reactants <= 2
