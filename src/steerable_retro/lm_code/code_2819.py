#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects a linear synthesis strategy where each reaction
    has only one non-commercial product that feeds into the next step.
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            # Count non-commercial molecule children
            non_commercial_children = 0
            for child in node.get("children", []):
                if child["type"] == "mol" and not child.get("in_stock", False):
                    non_commercial_children += 1

            # If more than one non-commercial child, it's not linear
            if non_commercial_children > 1:
                is_linear = False
                print(
                    f"Non-linear synthesis detected at depth {depth}: {non_commercial_children} non-commercial children"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Must have at least 3 reactions to be considered a meaningful linear strategy
    if reaction_count < 3:
        is_linear = False
        print(
            f"Only {reaction_count} reactions found, not enough for linear strategy classification"
        )

    return is_linear
