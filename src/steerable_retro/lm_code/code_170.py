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
    Detects if the synthetic route follows a linear strategy (vs convergent).
    Checks if most reactions have only 1-2 reactants.
    """
    reaction_count = 0
    linear_reaction_count = 0

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reaction_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]

            # Count number of reactants
            reactants = reactants_str.split(".")
            reaction_count += 1

            # If reaction has 1-2 reactants, consider it linear
            if len(reactants) <= 2:
                linear_reaction_count += 1

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If more than 80% of reactions are linear, consider the whole synthesis linear
    if reaction_count > 0:
        linear_ratio = linear_reaction_count / reaction_count
        print(f"Linear reactions: {linear_reaction_count}/{reaction_count} ({linear_ratio:.2f})")
        return linear_ratio >= 0.8

    return False
