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
    This function detects if the synthesis follows a linear strategy where each reaction
    has exactly one non-commercial reactant.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "children" in node:
                # Count non-commercial reactants (those that have children)
                non_commercial_count = sum(
                    1
                    for child in node["children"]
                    if child["type"] == "mol" and not child.get("in_stock", False)
                )

                if non_commercial_count > 1:
                    is_linear = False
                    print(
                        f"Found non-linear step at depth {depth} with {non_commercial_count} non-commercial reactants"
                    )

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if is_linear:
        print("Synthesis follows a linear strategy")

    return is_linear
