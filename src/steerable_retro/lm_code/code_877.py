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
    Detects if the synthesis follows a linear strategy where each reaction
    builds upon the product of the previous reaction.
    """
    reaction_count = 0
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, is_linear

        if node["type"] == "reaction":
            reaction_count += 1

            # Check if this reaction has exactly one non-commercial reactant
            # (which would be the product of the previous reaction)
            non_commercial_reactants = 0
            for child in node.get("children", []):
                if child["type"] == "mol":
                    if not child.get("in_stock", False):
                        non_commercial_reactants += 1

            # If more than one non-commercial reactant, it's not a linear synthesis
            if non_commercial_reactants > 1:
                is_linear = False
                print(
                    f"Non-linear step detected at depth {depth}: {non_commercial_reactants} non-commercial reactants"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # A synthesis with at least 3 reactions that maintains linearity
    if reaction_count >= 3 and is_linear:
        print(f"Linear synthesis strategy detected with {reaction_count} reactions")
        return True
    return False
