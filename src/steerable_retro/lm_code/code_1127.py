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
    This function detects if the synthesis follows a linear (non-convergent) strategy.

    A linear synthesis is characterized by:
    1. Each reaction step has at most one non-in-stock reactant
    2. There is a single path from the target molecule to the starting materials
    """
    # Track if we have a linear synthesis
    is_linear = True

    def count_non_stock_reactants(node):
        """Count the number of non-in-stock reactants in a reaction node"""
        non_stock_count = 0
        for child in node.get("children", []):
            if child["type"] == "mol" and not child.get("in_stock", False):
                non_stock_count += 1
        return non_stock_count

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        # If we already determined it's not linear, no need to continue
        if not is_linear:
            return

        # For reaction nodes, check if there's more than one non-in-stock reactant
        if node["type"] == "reaction":
            non_stock_count = count_non_stock_reactants(node)
            if non_stock_count > 1:
                print(
                    f"Found non-linear step: reaction has {non_stock_count} non-in-stock reactants"
                )
                is_linear = False
                return

        # Continue traversing children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the target molecule
    dfs_traverse(route)

    print(f"Synthesis is {'linear' if is_linear else 'convergent'}")
    return is_linear
