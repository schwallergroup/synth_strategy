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
    Detects if the synthetic route contains multiple fragment coupling reactions (3 or more).
    """
    fragment_coupling_count = 0

    def dfs_traverse(node):
        nonlocal fragment_coupling_count

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count as fragment coupling if there are 2 or more reactants
                if len(reactants) >= 2:
                    fragment_coupling_count += 1
                    print(
                        f"Found fragment coupling reaction (count: {fragment_coupling_count})"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if there are 3 or more fragment couplings
    return fragment_coupling_count >= 3
