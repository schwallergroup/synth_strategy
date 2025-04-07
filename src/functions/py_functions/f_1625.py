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
    Detects if the synthesis follows a linear strategy without convergent steps.
    """
    is_linear = True
    max_children_per_reaction = 0

    def dfs_traverse(node):
        nonlocal is_linear, max_children_per_reaction

        if node["type"] == "reaction":
            # Count the number of reactant children
            reactant_count = sum(
                1 for child in node.get("children", []) if child["type"] == "mol"
            )
            max_children_per_reaction = max(max_children_per_reaction, reactant_count)

            # If any reaction has more than 2 reactants, it's potentially convergent
            if reactant_count > 2:
                is_linear = False
                print(
                    f"Found potential convergent step with {reactant_count} reactants"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # A truly linear synthesis typically has at most 2 reactants per step
    if is_linear and max_children_per_reaction <= 2:
        print("Synthesis follows a linear strategy")
        return True
    return False
