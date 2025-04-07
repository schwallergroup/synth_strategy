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
    This function detects if the synthesis follows a linear strategy without convergent steps.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Count molecule children (reactants)
            mol_children = [
                child for child in node.get("children", []) if child["type"] == "mol"
            ]

            # If more than one significant reactant (not just a reagent), it's convergent
            significant_reactants = [
                child for child in mol_children if not child.get("in_stock", False)
            ]
            if len(significant_reactants) > 1:
                is_linear = False
                print(
                    f"Convergent step detected with {len(significant_reactants)} significant reactants"
                )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    if is_linear:
        print("Linear synthesis strategy confirmed")

    return is_linear
