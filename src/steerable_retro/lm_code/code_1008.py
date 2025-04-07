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
    Detects if the synthesis follows a linear strategy where each reaction
    has only one complex product that is used in the next step.
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Check if this reaction has multiple complex products
            complex_product_count = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    complex_product_count += 1

            if complex_product_count > 1:
                is_linear = False
                print(f"Non-linear step detected at depth {node.get('depth', 'unknown')}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # A synthesis with only one reaction is trivially linear
    if reaction_count <= 1:
        is_linear = False

    return is_linear and reaction_count > 1
