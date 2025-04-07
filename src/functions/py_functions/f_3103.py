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
    Detects if the route follows a linear fragment assembly strategy where
    each step adds one fragment at a time rather than converging multiple fragments.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If any reaction has more than 2 reactants, it's not a simple linear assembly
            if len(reactants) > 2:
                print(
                    f"Found reaction with {len(reactants)} reactants, not linear assembly"
                )
                is_linear = False

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if is_linear:
        print("Route follows linear fragment assembly strategy")

    return is_linear
