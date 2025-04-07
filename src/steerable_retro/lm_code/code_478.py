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
    This function detects if the synthesis follows a linear strategy (no convergent steps).
    A linear synthesis has exactly one non-starting material reactant in each step.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Get the children of this reaction node (the reactants)
            reactant_nodes = node.get("children", [])

            # Count non-starting material reactants
            non_starting_material_count = 0
            for child in reactant_nodes:
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_starting_material_count += 1

            # If more than one non-starting material reactant, it's not a linear synthesis
            if non_starting_material_count > 1:
                is_linear = False
                print(
                    f"Non-linear step detected with {non_starting_material_count} non-starting material reactants"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Linear synthesis strategy: {is_linear}")
    return is_linear
