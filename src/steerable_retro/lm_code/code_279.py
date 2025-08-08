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
    This function detects a linear synthesis strategy (as opposed to convergent)
    where each reaction has only one product that feeds into the next reaction.

    In retrosynthetic analysis, this means each molecule node (except possibly
    starting materials) should have at most one reaction child.
    """
    # Track the linearity of the synthesis
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count

        if node["type"] == "mol" and not node.get("in_stock", False):
            # For molecule nodes that aren't starting materials, check if they have
            # exactly one reaction child (linear synthesis)
            reaction_children = [
                child for child in node.get("children", []) if child["type"] == "reaction"
            ]

            if len(reaction_children) > 1:
                is_linear = False
                print(
                    f"Found non-linear branch: molecule has {len(reaction_children)} reaction paths"
                )

            if len(reaction_children) == 1:
                # Only count reactions that are part of the main synthetic path
                reaction_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # A synthesis must have at least 2 reactions to be considered linear
    strategy_present = is_linear and reaction_count >= 2

    if strategy_present:
        print("Linear synthesis strategy detected")
    else:
        print("Linear synthesis strategy not detected")
        print(f"Is linear: {is_linear}, Reaction count: {reaction_count}")

    return strategy_present
