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
    Detects if the synthesis route involves a convergent strategy with multiple fragment couplings.
    """
    fragment_couplings = 0

    def dfs_traverse(node):
        nonlocal fragment_couplings

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of distinct reactant fragments
                if len(reactants) >= 2:
                    fragment_couplings += 1
                    print(f"Fragment coupling detected with {len(reactants)} reactants: {rsmi}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    result = fragment_couplings >= 2  # At least 2 fragment couplings for convergent synthesis
    print(
        f"Convergent fragment coupling strategy detected: {result} (found {fragment_couplings} couplings)"
    )
    return result
