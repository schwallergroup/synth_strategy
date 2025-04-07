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
    Detects if the synthesis follows a linear strategy (as opposed to convergent).
    Linear synthesis typically has one main intermediate that is sequentially modified.
    """
    # Track the number of reactants in each step
    reactant_counts = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                reactant_counts.append((len(reactants), depth))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort by depth
    reactant_counts.sort(key=lambda x: x[1])

    # Check if most steps have exactly 2 reactants (typical for linear synthesis)
    two_reactant_steps = sum(1 for count, _ in reactant_counts if count == 2)

    # If more than 75% of steps have 2 reactants, consider it linear
    is_linear = two_reactant_steps >= 0.75 * len(reactant_counts)

    if is_linear:
        print("Detected linear synthesis strategy")

    return is_linear
