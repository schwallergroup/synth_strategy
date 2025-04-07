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
    This function detects if the synthesis follows a linear strategy rather than a convergent one.
    A linear synthesis has a single branch at each step.
    """
    is_linear = True
    max_branches = 1

    def count_branches(node):
        nonlocal is_linear, max_branches

        if node["type"] == "reaction":
            # Count the number of reactant branches
            reactant_count = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    reactant_count += 1

            if reactant_count > max_branches:
                max_branches = reactant_count

            if reactant_count > 1:
                is_linear = False
                print(f"Found convergent step with {reactant_count} non-stock reactants")

        # Traverse children
        for child in node.get("children", []):
            count_branches(child)

    # Start traversal
    count_branches(route)

    print(f"Linear synthesis strategy detected: {is_linear} (max branches: {max_branches})")
    return is_linear
