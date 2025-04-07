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
    Detects if the synthesis follows a linear strategy (as opposed to convergent).
    Linear synthesis typically has reactions with 1-2 reactants throughout.
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count non-empty reactants
            reactant_count = sum(1 for r in reactants if r.strip())

            # If more than 2 reactants, likely not linear
            if reactant_count > 2:
                is_linear = False
                print(f"Found non-linear step with {reactant_count} reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Must have at least 3 reactions to be considered a meaningful linear synthesis
    strategy_present = is_linear and reaction_count >= 3
    print(
        f"Detected {reaction_count} reactions in a {'linear' if is_linear else 'non-linear'} sequence"
    )
    print(f"Strategy detection result: {strategy_present}")

    return strategy_present
