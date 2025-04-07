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
    This function detects if the synthesis follows a linear strategy (no convergent steps).

    A linear synthesis has each reaction step with only one non-in-stock reactant.
    A convergent synthesis has at least one step with multiple non-in-stock reactants coming together.

    Args:
        route: A synthesis route following the SynthesisRoute JSON schema

    Returns:
        bool: True if the synthesis is linear, False if it's convergent
    """
    # Track the maximum branching factor of non-in-stock reactants
    max_non_stock_children = 0

    def dfs_traverse(node):
        nonlocal max_non_stock_children

        if node["type"] == "reaction":
            # Count number of non-in-stock reactant children
            non_stock_reactants = sum(
                1
                for child in node.get("children", [])
                if child["type"] == "mol" and not child.get("in_stock", False)
            )

            max_non_stock_children = max(max_non_stock_children, non_stock_reactants)

            print(
                f"Reaction node with {len(node.get('children', []))} total children and {non_stock_reactants} non-in-stock reactants"
            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If max_non_stock_children is 1 or 0, it's a linear synthesis
    # If max_non_stock_children > 1, there's at least one convergent step
    is_linear = max_non_stock_children <= 1

    # Debug information
    print(f"Maximum non-in-stock reactants in any step: {max_non_stock_children}")
    if is_linear:
        print("Detected linear synthesis strategy (no convergent steps)")
    else:
        print(
            f"Detected convergent synthesis with maximum {max_non_stock_children} non-in-stock reactants in a step"
        )

    return is_linear
