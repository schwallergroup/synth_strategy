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
    This function detects if the synthesis follows a linear strategy with no branching.

    In a linear synthesis:
    - Each molecule node should lead to at most one reaction node
    - Reaction nodes can have multiple children (reactants) without affecting linearity
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        # Only check branching for molecule nodes
        if node.get("type") == "mol" and not node.get("in_stock", False):
            reaction_children = [
                child for child in node.get("children", []) if child.get("type") == "reaction"
            ]

            # If a molecule leads to more than one reaction, it's not linear
            if len(reaction_children) > 1:
                print(f"Found branching at molecule: {node.get('smiles', 'unknown')}")
                is_linear = False

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Linear synthesis strategy: {is_linear}")
    return is_linear
