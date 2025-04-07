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
    Detects if the synthesis follows a linear heterocycle elaboration strategy
    with sequential functionalization rather than convergent synthesis.
    """
    # Track reactions and their depths
    reactions_by_depth = {}

    def dfs_traverse(node):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            depth = node.get("depth", -1)
            if depth not in reactions_by_depth:
                reactions_by_depth[depth] = []
            reactions_by_depth[depth].append(node)

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a linear sequence (one reaction per depth)
    depths = sorted(reactions_by_depth.keys())
    is_linear = all(len(reactions_by_depth[d]) == 1 for d in depths)

    # Check if we have at least 4 steps
    if is_linear and len(depths) >= 4:
        print("Detected linear heterocycle elaboration strategy")
        return True
    return False
