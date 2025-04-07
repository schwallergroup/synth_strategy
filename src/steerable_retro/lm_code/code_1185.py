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
    This function detects if the synthetic route follows a linear synthesis pattern
    rather than a convergent one.
    """
    # Track the maximum branching factor in the synthesis tree
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            # Count the number of reactants in this reaction
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                num_reactants = len([r for r in reactants if r])
                max_branching = max(max_branching, num_reactants)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If max branching is <= 2, it's likely a linear synthesis
    # (allowing for one main reactant and one reagent)
    is_linear = max_branching <= 2
    print(f"Maximum branching factor: {max_branching}, Linear synthesis: {is_linear}")
    return is_linear
