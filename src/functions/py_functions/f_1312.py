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
    Detects a linear synthesis strategy (as opposed to convergent) by checking if most reactions
    have only 1-2 reactants.
    """
    reaction_count = 0
    linear_reaction_count = 0

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reaction_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            reaction_count += 1
            if len(reactants) <= 2:
                linear_reaction_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If at least 80% of reactions are linear (1-2 reactants), consider it a linear strategy
    if reaction_count > 0 and (linear_reaction_count / reaction_count) >= 0.8:
        print(
            f"Linear synthesis strategy detected: {linear_reaction_count}/{reaction_count} reactions are linear"
        )
        return True
    return False
