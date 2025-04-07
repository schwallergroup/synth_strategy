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
    Detects if the synthesis route follows a linear fragment assembly approach
    without convergent steps.
    """
    # Track the maximum number of reactants in any step
    max_reactants = 0

    def dfs_traverse(node):
        nonlocal max_reactants

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count number of reactants
            num_reactants = len(reactants)
            max_reactants = max(max_reactants, num_reactants)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If max_reactants is consistently 2 or less, it's likely a linear synthesis
    is_linear = max_reactants <= 2
    print(f"Maximum reactants in any step: {max_reactants}, Linear synthesis: {is_linear}")
    return is_linear
