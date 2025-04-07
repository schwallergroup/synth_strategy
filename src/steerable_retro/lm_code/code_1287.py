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
    This function detects a linear fragment assembly strategy where fragments are added sequentially.
    It checks for a pattern of sequential reactions where each step adds one fragment.
    """
    # Track the number of fragments joined
    fragments_joined = 0
    linear_strategy = False

    def dfs_traverse(node):
        nonlocal fragments_joined, linear_strategy

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If reaction has multiple reactants, it's joining fragments
            if len(reactants) > 1:
                fragments_joined += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If we have a sequence of fragment-joining reactions (at least 2)
    # and they occur in a linear fashion, it's a linear fragment assembly
    if fragments_joined >= 2:
        print(f"Linear fragment assembly detected with {fragments_joined} fragment-joining steps")
        linear_strategy = True

    return linear_strategy
