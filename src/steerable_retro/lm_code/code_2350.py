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
    This function detects a linear synthesis strategy as opposed to a convergent one.
    It checks if most reactions have only one non-commercial reactant.
    """
    reaction_count = 0
    linear_reaction_count = 0

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Get children nodes (reactants)
            children = node.get("children", [])

            # Count non-commercial reactants
            non_commercial_reactants = 0
            for child in children:
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_commercial_reactants += 1

            # If only one non-commercial reactant, it's a linear step
            if non_commercial_reactants <= 1:
                linear_reaction_count += 1
                print(f"Found linear reaction step")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If more than 70% of reactions are linear, consider it a linear synthesis
    if reaction_count > 0 and linear_reaction_count / reaction_count >= 0.7:
        return True
    return False
