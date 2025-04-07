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
    This function detects if the synthetic route follows a linear synthesis approach
    (each step has only one complex reactant from previous step).
    """
    # For linear synthesis, we need to check that at each reaction step,
    # there's only one non-commercial reactant

    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            complex_reactants = 0

            # Count children that are molecules and not in_stock
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    complex_reactants += 1

            # If more than one complex reactant, it's not linear
            if complex_reactants > 1:
                print(
                    f"Non-linear step detected with {complex_reactants} complex reactants"
                )
                is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return is_linear
