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
    Detects a linear fragment assembly strategy where each step adds one fragment
    without convergent synthesis patterns.
    """
    # Track the number of reactions and whether any have more than 2 reactants
    reaction_count = 0
    has_convergent_step = False

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, has_convergent_step

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = [r for r in rsmi.split(">")[0].split(".") if r]

            # If more than 2 reactants, might be convergent synthesis
            if len(reactants) > 2:
                has_convergent_step = True
                print(
                    f"Found potential convergent step at depth {depth} with {len(reactants)} reactants"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Linear strategy has multiple reactions but no convergent steps
    strategy_present = reaction_count >= 3 and not has_convergent_step

    print(
        f"Strategy detection result: {strategy_present} (reactions: {reaction_count}, convergent: {has_convergent_step})"
    )
    return strategy_present
