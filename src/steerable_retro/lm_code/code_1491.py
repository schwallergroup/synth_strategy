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
    Detects if the synthesis follows a linear pattern where each step builds directly
    on the previous product.
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Check if this reaction has exactly one product
            rsmi = node["metadata"]["rsmi"]
            products = rsmi.split(">")[-1].split(".")

            if len(products) > 1:
                is_linear = False
                print("Non-linear synthesis detected: reaction has multiple products")

            # Check if this reaction has more than 2 reactants (suggesting convergent synthesis)
            reactants = rsmi.split(">")[0].split(".")
            if len(reactants) > 2:
                # This is a simplification - some linear syntheses might use multiple reagents
                # A more sophisticated check would distinguish between reagents and key building blocks
                is_linear = False
                print("Potentially convergent synthesis detected: reaction has multiple reactants")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Need at least 2 reactions to determine if synthesis is linear
    if reaction_count < 2:
        return False

    return is_linear
