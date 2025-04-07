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
    This function detects a linear synthesis strategy (as opposed to convergent)
    by checking if each reaction has only one non-commercial reactant.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count the number of non-commercial reactants
            non_commercial_reactants = 0

            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_commercial_reactants += 1

            # If more than one non-commercial reactant, it's not a linear synthesis
            if non_commercial_reactants > 1:
                print(
                    f"Found convergent step at depth {depth} with {non_commercial_reactants} non-commercial reactants"
                )
                is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if is_linear:
        print("Confirmed linear synthesis strategy")

    return is_linear
