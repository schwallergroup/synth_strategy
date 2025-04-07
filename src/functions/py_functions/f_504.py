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
    This function detects if the synthetic route follows a linear strategy (as opposed to convergent).
    Linear synthesis typically has 1-2 reactants per step.
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            reaction_count += 1

            # If more than 2 reactants, it might be a convergent step
            if len(reactants) > 2:
                is_linear = False
                print(f"Convergent step detected with {len(reactants)} reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Only consider routes with multiple reactions
    if reaction_count >= 3:
        print(
            f"Route has {reaction_count} reactions and is {'linear' if is_linear else 'convergent'}"
        )
        return is_linear

    return False
