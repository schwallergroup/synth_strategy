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
    Detects if the route follows a predominantly linear synthesis strategy
    with sequential transformations until a late-stage coupling.
    """
    # Track reactions with their depth
    reactions_with_depth = []

    def count_reactants(node):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            return len(reactants)
        return 0

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            branching_factor = count_reactants(node)
            reactions_with_depth.append((depth, branching_factor))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (lowest depth = late-stage, highest depth = early-stage)
    reactions_with_depth.sort(key=lambda x: x[0])

    # Extract just the branching factors in order from late-stage to early-stage
    branching_factors = [bf for _, bf in reactions_with_depth]

    print(f"Reactions with depth: {reactions_with_depth}")
    print(f"Branching factors (late to early): {branching_factors}")

    # Check if we have enough reactions to analyze
    if len(branching_factors) < 2:
        print("Not enough reactions to determine strategy")
        return False

    # Check if the late-stage reaction (first in our sorted list) has multiple reactants
    late_stage_coupling = branching_factors[0] > 1

    # Check if most of the early-stage reactions are linear (one reactant)
    early_reactions = branching_factors[1:]
    linear_count = sum(1 for bf in early_reactions if bf == 1)
    mostly_linear = linear_count >= len(early_reactions) * 0.7

    is_linear_with_late_coupling = late_stage_coupling and mostly_linear

    if is_linear_with_late_coupling:
        print("Found linear synthesis with late-stage coupling")
        print(f"Late-stage coupling has {branching_factors[0]} reactants")
        print(f"{linear_count} out of {len(early_reactions)} early-stage reactions are linear")
    else:
        if not late_stage_coupling:
            print("Late-stage reaction is not a coupling (only one reactant)")
        if not mostly_linear:
            print("Early-stage reactions are not predominantly linear")

    return is_linear_with_late_coupling
