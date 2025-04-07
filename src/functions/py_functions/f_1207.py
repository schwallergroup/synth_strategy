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
    Detects if the route follows a linear synthesis strategy with no convergent steps.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count the number of non-stock reactants
            non_stock_reactants = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_stock_reactants += 1

            # If more than one non-stock reactant, it's a convergent step
            if non_stock_reactants > 1:
                is_linear = False
                print(
                    f"Found convergent step with {non_stock_reactants} non-stock reactants"
                )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Linear synthesis strategy detected: {is_linear}")
    return is_linear
