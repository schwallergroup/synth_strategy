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
    Detects if the synthesis follows a linear strategy (no convergent steps).
    """
    # In a linear synthesis, each reaction node should have exactly one molecule child
    # that is not a starting material

    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            non_starting_material_children = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_starting_material_children += 1

            if non_starting_material_children > 1:
                is_linear = False
                print(
                    f"Convergent step detected: {non_starting_material_children} non-starting material reactants"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if is_linear:
        print("Linear synthesis strategy confirmed")

    return is_linear
