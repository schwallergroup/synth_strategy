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
    This function detects if the synthesis follows a linear strategy
    (as opposed to convergent) by checking if most reactions have only 1-2 reactants.
    """
    reaction_count = 0
    linear_reaction_count = 0

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reaction_count

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if rsmi:
                reaction_count += 1
                reactants = rsmi.split(">")[0].split(".")
                # Count non-empty reactants
                reactant_count = sum(1 for r in reactants if r.strip())

                # Linear reactions typically have 1-2 reactants
                if reactant_count <= 2:
                    linear_reaction_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If at least 80% of reactions are linear, consider it a linear synthesis
    if reaction_count > 0 and (linear_reaction_count / reaction_count) >= 0.8:
        print(
            f"Linear synthesis detected: {linear_reaction_count}/{reaction_count} reactions are linear"
        )
        return True
    return False
