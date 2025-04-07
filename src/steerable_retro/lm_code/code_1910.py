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
    Detects if the synthetic route follows a linear strategy without convergent steps.

    A linear synthesis is one where each reaction has at most one non-starting material
    reactant. This means the synthesis proceeds in a straight line without convergent steps.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and is_linear:
            # Count non-starting material reactants
            non_starting_material_count = 0

            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_starting_material_count += 1

            # If more than one non-starting material, it's not linear
            if non_starting_material_count > 1:
                is_linear = False
                print(
                    f"Found non-linear reaction with {non_starting_material_count} non-starting material reactants"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Is linear synthesis: {is_linear}")

    return is_linear
