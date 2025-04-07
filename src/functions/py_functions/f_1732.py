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
    This function detects if the synthesis follows a linear strategy (each step builds on the previous
    without convergent steps).
    """
    is_linear = True

    def count_children(node):
        if node["type"] == "reaction":
            children = node.get("children", [])
            # If a reaction has more than 2 children (more than 1 non-trivial reactant),
            # it's likely a convergent step
            non_trivial_children = 0
            for child in children:
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_trivial_children += 1
                elif child["type"] == "reaction":
                    non_trivial_children += 1

            if non_trivial_children > 1:
                return False

            # Recursively check children
            for child in children:
                if not count_children(child):
                    return False

        return True

    is_linear = count_children(route)
    if is_linear:
        print("Linear synthesis strategy detected")
    else:
        print("Convergent synthesis strategy detected")

    return is_linear
