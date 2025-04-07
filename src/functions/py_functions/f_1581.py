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
    Detects if the synthetic route follows a linear build-up strategy without convergent steps.
    """
    # In a linear synthesis, each reaction node should have exactly one molecule child
    # that is not a starting material
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count non-starting material children
            non_starting_material_children = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_starting_material_children += 1

            # If more than one non-starting material child, it's not linear
            if non_starting_material_children > 1:
                is_linear = False
                print("Detected non-linear (convergent) synthesis")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if is_linear:
        print("Detected linear synthesis strategy")

    return is_linear
