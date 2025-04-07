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
    This function detects a linear synthesis strategy with sequential fragment additions
    rather than convergent synthesis.
    """
    max_fragments_per_reaction = 0

    def dfs_traverse(node):
        nonlocal max_fragments_per_reaction

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count number of fragments being combined
                num_fragments = len(reactants)
                max_fragments_per_reaction = max(
                    max_fragments_per_reaction, num_fragments
                )

                if num_fragments > 2:
                    print(f"Found reaction with {num_fragments} fragments")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # If max fragments per reaction is 2, it's likely a linear synthesis
    return max_fragments_per_reaction <= 2
