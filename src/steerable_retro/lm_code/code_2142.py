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
    This function detects if the synthesis follows a linear strategy (each step typically has 1-2 reactants).
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If more than 2 reactants, likely not linear
                if len(reactants) > 2:
                    is_linear = False
                    print(f"Non-linear step detected with {len(reactants)} reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Must have at least 3 reactions to be considered a meaningful linear synthesis
    return is_linear and reaction_count >= 3
