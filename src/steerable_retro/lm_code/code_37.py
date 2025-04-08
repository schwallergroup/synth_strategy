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
    Detects if the synthesis uses a strategy of combining multiple fragments (3 or more)
    throughout the synthesis route.
    """
    fragment_combination_count = 0

    def dfs_traverse(node):
        nonlocal fragment_combination_count

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_part = rsmi.split(">")[0]
            reactants = reactants_part.split(".")

            # If we have 2 or more reactants, this is a fragment combination
            if len(reactants) >= 2:
                fragment_combination_count += 1
                print(
                    f"Detected fragment combination at depth {node.get('metadata', {}).get('depth', -1)}"
                )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Return True if we have 3 or more fragment combinations
    return fragment_combination_count >= 3
