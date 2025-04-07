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
    This function detects if the synthesis follows a linear strategy
    without convergent steps (each reaction has only one product molecule).
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains multiple molecules (indicated by ".")
                if "." in product_smiles:
                    # This might be a convergent step or a reaction with byproducts
                    # For simplicity, we'll consider it non-linear
                    is_linear = False
                    print(f"Non-linear step detected at depth {depth}: multiple products")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if is_linear:
        print("Linear synthesis strategy confirmed")
    return is_linear
