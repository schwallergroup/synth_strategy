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
    This function detects if the synthetic route follows a linear synthesis strategy
    (each reaction builds upon a single previous product).
    """
    is_linear = True
    max_branching = 0

    def count_mol_children(node):
        if node["type"] == "reaction":
            mol_children = sum(
                1 for child in node.get("children", []) if child["type"] == "mol"
            )
            return mol_children
        return 0

    def dfs_traverse(node):
        nonlocal is_linear, max_branching

        if node["type"] == "reaction":
            branching = count_mol_children(node)
            max_branching = max(max_branching, branching)

            if branching > 2:  # More than 2 reactants indicates convergent synthesis
                is_linear = False
                print(f"Found convergent step with {branching} reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Linear synthesis strategy detected: {is_linear} (max branching: {max_branching})"
    )
    return is_linear
