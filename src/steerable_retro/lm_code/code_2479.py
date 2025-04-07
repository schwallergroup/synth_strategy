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
    This function detects if the synthesis follows a linear fragment assembly strategy
    with at least 4 distinct fragments being added sequentially.
    """
    fragment_count = 0
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal fragment_count, max_depth

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If there are multiple reactants, it's potentially a fragment addition
                if len(reactants) > 1:
                    fragment_count += 1
                    print(f"Fragment addition detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have at least 4 fragment additions and the synthesis is mostly linear
    is_linear_assembly = fragment_count >= 4
    print(
        f"Linear fragment assembly detected: {is_linear_assembly} with {fragment_count} fragments"
    )
    return is_linear_assembly
