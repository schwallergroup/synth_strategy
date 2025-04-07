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
    This function detects if the synthesis follows a linear rather than convergent strategy.
    """
    max_branch_factor = 0

    def dfs_traverse(node):
        nonlocal max_branch_factor

        if node["type"] == "reaction":
            # Count number of reactants in this reaction
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                branch_factor = len(reactants)
                max_branch_factor = max(max_branch_factor, branch_factor)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If max branch factor is 1 or 2, consider it a linear synthesis
    # (allowing for 2 to account for reagents that aren't part of the main skeleton)
    is_linear = max_branch_factor <= 2
    print(f"Maximum branch factor: {max_branch_factor}, Linear synthesis: {is_linear}")
    return is_linear
