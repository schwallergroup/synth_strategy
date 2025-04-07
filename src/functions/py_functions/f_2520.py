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
    This function detects if the synthetic route follows a linear strategy without convergent steps.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count non-commercial reactants
            non_commercial_reactants = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_commercial_reactants += 1

            # If more than one non-commercial reactant, it's convergent
            if non_commercial_reactants > 1:
                is_linear = False
                print(
                    "Convergent step detected with",
                    non_commercial_reactants,
                    "non-commercial reactants",
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Linear synthesis strategy detected: {is_linear}")
    return is_linear
