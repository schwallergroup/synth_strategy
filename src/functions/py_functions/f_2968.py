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
    Detects if the synthesis follows a linear strategy with a late-stage coupling of fragments.
    """
    late_coupling = False
    linear_steps = 0
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal late_coupling, linear_steps, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Check if this is a coupling reaction (combining multiple fragments)
                if len(reactants) >= 2 and depth <= 1:  # Late stage
                    # Look for indicators of coupling reactions
                    if "#" in rsmi:  # Alkyne coupling
                        late_coupling = True
                        print(f"Late-stage coupling detected at depth {depth}")

            linear_steps += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # A synthesis is considered linear if it has multiple steps and ends with a coupling
    return late_coupling and linear_steps >= 3 and max_depth >= 3
