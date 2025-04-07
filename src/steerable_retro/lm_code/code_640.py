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
    This function detects a linear fragment assembly strategy where the molecule
    is built up sequentially without convergent steps.
    """
    # Track branching factor at each node
    branching_factors = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            # Count number of reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            num_reactants = len(reactants)
            branching_factors.append(num_reactants)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if all reactions have at most 2 reactants (linear assembly)
    is_linear = all(bf <= 2 for bf in branching_factors)
    print(f"Branching factors: {branching_factors}, Linear: {is_linear}")

    return (
        is_linear and len(branching_factors) >= 3
    )  # At least 3 reactions for a meaningful strategy
