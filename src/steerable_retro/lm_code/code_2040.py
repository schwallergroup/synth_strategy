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
    This function detects if the synthesis follows a predominantly linear path
    rather than a convergent approach.
    """
    # Track the maximum branching factor in the synthesis tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count the number of reactants in this reaction
            children = node.get("children", [])
            reactant_count = sum(1 for child in children if child["type"] == "mol")
            max_branching = max(max_branching, reactant_count)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If the maximum branching factor is <= 2, consider it a linear synthesis
    result = max_branching <= 2
    if result:
        print("Linear synthesis strategy detected")
    return result
