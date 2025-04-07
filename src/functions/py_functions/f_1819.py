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
    Detects if the synthesis follows a linear pattern (no convergent steps).

    A linear synthesis has no convergent steps, meaning each reaction has at most
    one non-stock (synthesized) reactant. In the retrosynthetic tree, this means
    each reaction node has at most one child that is a non-stock molecule.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count non-stock molecule children
            non_stock_children = 0

            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_stock_children += 1

            # If more than one non-stock child, it's a convergent synthesis
            if non_stock_children > 1:
                is_linear = False
                print(
                    f"Found convergent step with {non_stock_children} non-stock reactants"
                )

        # Continue traversal if still potentially linear
        if is_linear:
            for child in node.get("children", []):
                dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return is_linear
