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
    Checks if the synthesis route is convergent (has multiple branches).
    A convergent synthesis has at least one reaction with multiple reactants.
    """
    # Track if we've found a convergent step
    found_convergent = False

    def dfs(node, depth=0):
        nonlocal found_convergent

        # Skip leaf nodes (starting materials)
        if node.get("type") == "mol" and node.get("in_stock", False):
            return

        # Check reaction nodes for convergent steps
        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            # If there are multiple reactants (separated by "."), it's convergent
            if "." in reactants_part:
                print(f"Found convergent step: {rsmi}")
                found_convergent = True

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS from the root
    dfs(route)
    return found_convergent
