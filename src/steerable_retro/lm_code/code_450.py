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
    This function detects if the synthetic route follows a linear strategy
    without convergent steps.
    """
    is_linear = True

    def count_reactants(reaction_node):
        rsmi = reaction_node["metadata"].get("rsmi", "")
        if not rsmi:
            return 0

        reactants_part = rsmi.split(">")[0]
        # Count number of distinct reactants (separated by ".")
        return len(reactants_part.split("."))

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # If any reaction has more than 2 reactants, it might be convergent
            # (allowing for 2 because many reactions have a reagent in addition to the main substrate)
            if count_reactants(node) > 2:
                is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if is_linear:
        print("Linear synthesis strategy detected")
    return is_linear
