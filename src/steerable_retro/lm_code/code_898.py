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
    This function detects a linear synthetic strategy where complexity is built
    through sequential addition of fragments rather than convergent synthesis.
    """
    # Track reaction depths and branching
    reaction_depths = []
    max_children_per_node = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_children_per_node

        if node["type"] == "reaction":
            reaction_depths.append(depth)

            # Count number of reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            num_reactants = len(reactants)
            max_children_per_node = max(max_children_per_node, len(node.get("children", [])))

            print(
                f"Reaction at depth {depth} has {num_reactants} reactants and {len(node.get('children', []))} children"
            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze reaction pattern
    # Linear synthesis typically has:
    # 1. Reactions at multiple depths (sequential)
    # 2. Limited branching (max 2 children per node)
    # 3. At least 3 reaction steps

    is_linear = (
        len(reaction_depths) >= 3
        and max_children_per_node <= 2  # At least 3 reaction steps
        and len(set(reaction_depths)) >= 3  # Limited branching  # Reactions at multiple depths
    )

    print(f"Strategy detection results:")
    print(f"  - Number of reaction steps: {len(reaction_depths)}")
    print(f"  - Max children per node: {max_children_per_node}")
    print(f"  - Unique reaction depths: {len(set(reaction_depths))}")
    print(f"  - Linear strategy detected: {is_linear}")

    return is_linear
