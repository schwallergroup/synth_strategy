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
    This function detects a linear synthesis strategy where each step
    adds one fragment at a time rather than converging multiple fragments.
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If more than 2 reactants, it's not a simple linear addition
                if len(reactants) > 2:
                    print(f"More than 2 reactants found at depth {depth}, not a linear strategy")
                    is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_linear
