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
    This function detects if the synthesis follows a linear build-up strategy
    rather than a convergent approach.

    A linear build-up strategy typically involves:
    1. A chain-like synthesis path with limited branching
    2. Sequential addition of building blocks
    3. Most reaction steps having 1-2 reactants
    """
    # Track metrics for analysis
    reaction_count = 0
    max_reactants_per_step = 0

    # Track tree structure
    max_depth = 0
    branch_points = 0

    def dfs_traverse(node, depth=0):
        nonlocal reaction_count, max_reactants_per_step, max_depth, branch_points

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of distinct reactants
                reactant_count = len([r for r in reactants if r])
                max_reactants_per_step = max(max_reactants_per_step, reactant_count)
                reaction_count += 1

                # Check for branching (more than 2 children indicates convergent synthesis)
                if len(node.get("children", [])) > 2:
                    branch_points += 1

        # Recursively process children
        children = node.get("children", [])
        for child in children:
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze synthesis pattern
    print(f"Reaction count: {reaction_count}")
    print(f"Max reactants per step: {max_reactants_per_step}")
    print(f"Max depth: {max_depth}")
    print(f"Branch points: {branch_points}")

    # Criteria for linear build-up:
    # 1. At least one reaction
    # 2. Limited branching (few branch points)
    # 3. Deeper tree structure (depth > branch points)
    if reaction_count > 0 and (branch_points == 0 or max_depth > branch_points * 2):
        print("Detected linear build-up strategy")
        return True

    return False
