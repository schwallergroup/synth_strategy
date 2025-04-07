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
    This function detects if the synthesis follows a linear build-up strategy
    without convergent steps.
    """
    max_reactants_per_step = 0

    def dfs_traverse(node):
        nonlocal max_reactants_per_step

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            num_reactants = len(reactants)
            max_reactants_per_step = max(max_reactants_per_step, num_reactants)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If max reactants per step is 2 or less, it's likely a linear build-up
    if max_reactants_per_step <= 2:
        print("Linear build-up strategy detected")
        return True
    return False
