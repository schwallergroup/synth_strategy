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
    This function detects a linear fragment assembly strategy as opposed to a convergent one.
    It checks if most reactions involve adding one fragment at a time rather than combining
    large fragments late in the synthesis.
    """
    reaction_count = 0
    linear_reaction_count = 0

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reaction_count

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            reaction_count += 1

            # If reaction has only one or two reactants, it's likely linear
            if len(reactants) <= 2:
                linear_reaction_count += 1
                print(f"Detected potential linear assembly step: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If most reactions are linear (one or two reactants), consider it a linear strategy
    return reaction_count > 0 and (linear_reaction_count / reaction_count) >= 0.7
