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
    Detects if the synthesis follows a linear build-up strategy rather than a convergent approach.
    Linear syntheses typically have fewer reactants per step.
    """
    is_linear = True
    steps_with_multiple_reactants = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, steps_with_multiple_reactants, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If a step has more than 2 reactants, it's less likely to be linear
            if len(reactants) > 2:
                steps_with_multiple_reactants += 1

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # If more than 25% of steps have multiple reactants, consider it less linear
    if total_steps > 0 and steps_with_multiple_reactants / total_steps > 0.25:
        is_linear = False

    print(f"Linear synthesis strategy detected: {is_linear}")
    print(f"Steps with multiple reactants: {steps_with_multiple_reactants} out of {total_steps}")

    return is_linear
