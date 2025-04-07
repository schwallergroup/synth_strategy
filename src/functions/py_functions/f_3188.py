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
    Detects if the synthesis follows a linear strategy (mostly one reactant per step).
    """
    linear_steps = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal linear_steps, total_steps

        if node["type"] == "reaction":
            total_steps += 1

            # Get reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If there are at most two reactants, consider it a linear step
            # This accounts for main reactant + simple reagent scenarios
            if reactants and len(reactants) <= 2:
                linear_steps += 1
                print(
                    f"Detected linear step {total_steps} with {len(reactants)} reactants: {reactants}"
                )
            else:
                print(
                    f"Non-linear step {total_steps} with {len(reactants)} reactants: {reactants}"
                )

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Linear steps: {linear_steps}, Total steps: {total_steps}")

    # Return True if most steps (>50%) are linear
    return total_steps > 0 and linear_steps / total_steps > 0.5
