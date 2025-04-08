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
    This function detects if the synthesis follows a linear strategy rather than
    a convergent one, by checking if most reactions have only 1-2 reactants.
    """
    reaction_count = 0
    linear_reaction_count = 0

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reaction_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                reaction_count += 1

                # If reaction has 1-2 reactants, consider it linear
                if len(reactants) <= 2:
                    linear_reaction_count += 1
                    print(f"Linear reaction detected: {len(reactants)} reactants")
                else:
                    print(f"Convergent reaction detected: {len(reactants)} reactants")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If at least 80% of reactions are linear, consider it a linear synthesis strategy
    if reaction_count > 0:
        linear_ratio = linear_reaction_count / reaction_count
        print(f"Linear ratio: {linear_ratio} ({linear_reaction_count}/{reaction_count})")
        return linear_ratio >= 0.8

    return False
