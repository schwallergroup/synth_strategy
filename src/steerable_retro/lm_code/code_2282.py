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
    This function detects if the synthesis employs a convergent strategy
    with multiple fragments being combined.
    """
    fragment_combination_count = 0

    def dfs_traverse(node):
        nonlocal fragment_combination_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If we have multiple reactants, it's a fragment combination
                if len(reactants) >= 2:
                    fragment_combination_count += 1
                    print(f"Found fragment combination: {rsmi}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if we found at least 2 fragment combinations
    result = fragment_combination_count >= 2
    print(f"Convergent synthesis strategy detected: {result} (count: {fragment_combination_count})")
    return result
