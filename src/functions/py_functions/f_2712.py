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
    Detects if the synthesis route is convergent (multiple fragments combined).
    """
    fragment_count = 0

    def dfs_traverse(node):
        nonlocal fragment_count

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If there are multiple reactants, this might be a fragment combination
            if len(reactants) > 1:
                # Check if reactants are complex enough to be considered fragments
                complex_reactants = [
                    r for r in reactants if r.count("c") > 3 or r.count("C") > 3
                ]
                if len(complex_reactants) > 1:
                    fragment_count += 1
                    print(f"Detected fragment combination: {rsmi}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return (
        fragment_count >= 2
    )  # At least 2 fragment combinations for convergent synthesis
