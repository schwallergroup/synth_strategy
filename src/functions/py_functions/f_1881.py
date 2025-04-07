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
    Detects if the synthesis follows a linear strategy rather than convergent.
    Linear synthesis is characterized by a chain-like reaction sequence.
    """
    # Count branching factor at each reaction node
    branching_factors = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            # Count number of reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            branching_factors.append(len(reactants))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If most reactions have only 1 or 2 reactants, it's likely a linear synthesis
    if (
        branching_factors
        and sum(1 for bf in branching_factors if bf <= 2) / len(branching_factors)
        >= 0.7
    ):
        print("Found linear synthesis strategy")
        return True
    return False
