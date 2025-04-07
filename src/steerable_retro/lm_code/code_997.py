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
    Detects if the synthesis follows a convergent approach where two or more fragments
    are combined in the early stages of the synthesis (high depth in retrosynthetic tree).
    """
    convergent_synthesis_found = False

    def dfs_traverse(node, depth=0):
        nonlocal convergent_synthesis_found

        if node["type"] == "reaction" and depth >= 3:  # Early in synthesis (high depth)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # If multiple reactants combine to form a single product
                if len(reactants) >= 2:
                    print(f"Convergent synthesis detected at depth {depth}")
                    convergent_synthesis_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return convergent_synthesis_found
