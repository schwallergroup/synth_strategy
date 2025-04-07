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
    Detects if the synthesis follows a linear strategy without convergent steps.
    Linear synthesis typically has exactly two reactants per step.
    """
    # Track if all steps are consistent with linear synthesis
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_count

        # Process reaction nodes
        if node["type"] == "reaction":
            reaction_count += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If more than 2 reactants, it might be convergent rather than linear
                if len(reactants) > 2:
                    print(
                        f"Found more than 2 reactants at depth {depth}, suggesting non-linear synthesis"
                    )
                    is_linear = False

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Ensure we found at least one reaction
    if reaction_count == 0:
        return False

    return is_linear
