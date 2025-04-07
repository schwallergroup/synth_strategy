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
    This function detects linear fragment assembly strategy (as opposed to convergent),
    where each reaction adds one fragment at a time.
    """
    linear_assembly = True

    def dfs_traverse(node):
        nonlocal linear_assembly

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If any reaction has more than 2 reactants, it's not a simple linear assembly
            if len(reactants) > 2:
                linear_assembly = False
                print("Found reaction with more than 2 reactants - not linear assembly")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return linear_assembly
