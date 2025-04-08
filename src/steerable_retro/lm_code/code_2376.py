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
    This function detects if the synthesis follows a convergent approach where
    multiple fragments are prepared separately and then joined.
    """
    fragment_count = 0
    joining_reactions = 0

    def dfs_traverse(node):
        nonlocal fragment_count, joining_reactions

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count distinct fragments in reactants
                if len(reactants) > 1:
                    fragment_count += len(reactants) - 1
                    joining_reactions += 1
                    print(f"Found joining reaction with {len(reactants)} fragments")

        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # A convergent synthesis typically has multiple fragments and joining reactions
    result = fragment_count >= 2 and joining_reactions >= 1
    print(
        f"Convergent synthesis strategy detected: {result} (fragments: {fragment_count}, joining reactions: {joining_reactions})"
    )
    return result
