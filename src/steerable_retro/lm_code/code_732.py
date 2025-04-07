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
    This function detects if the synthesis follows a linear strategy (as opposed to convergent).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count number of reactants that are not in_stock
            non_stock_reactants = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_stock_reactants += 1

            # If more than one non-stock reactant, it's not a linear synthesis
            if non_stock_reactants > 1:
                is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if is_linear:
        print("Detected linear synthesis strategy")

    return is_linear
