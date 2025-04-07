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
    This function detects if the synthesis follows a linear strategy without convergent steps.
    """
    is_linear = True
    max_children_per_reaction = 0

    def dfs_traverse(node):
        nonlocal is_linear, max_children_per_reaction

        if node["type"] == "reaction":
            # Count number of reactants (children)
            num_children = len(node.get("children", []))
            max_children_per_reaction = max(max_children_per_reaction, num_children)

            # If any reaction has more than 2 reactants, it might be convergent
            if num_children > 2:
                is_linear = False
                print(f"Found potential convergent step with {num_children} reactants")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # A truly linear synthesis typically has at most 2 reactants per step
    if max_children_per_reaction <= 2:
        print("Confirmed linear synthesis strategy")
    else:
        is_linear = False

    return is_linear
