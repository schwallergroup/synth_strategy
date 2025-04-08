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
    This function detects a convergent synthesis strategy where multiple fragments
    are combined in late-stage reactions.
    """
    has_convergent_step = False

    def dfs_traverse(node, depth=0):
        nonlocal has_convergent_step

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count non-empty reactants
            reactant_count = sum(1 for r in reactants if r.strip())

            # If this is a late-stage reaction (depth <= 1) with 3+ reactants
            if depth <= 1 and reactant_count >= 3:
                has_convergent_step = True
                print(
                    f"Detected convergent step at depth {depth} with {reactant_count} reactants: {rsmi}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Convergent synthesis strategy detected: {has_convergent_step}")
    return has_convergent_step
