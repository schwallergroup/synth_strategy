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
    This function detects linear synthesis pattern.
    Checks if the synthesis follows a linear path without convergent steps.
    """
    is_linear = True
    max_reactants_per_step = 2  # Allow up to 2 reactants per step for linear synthesis

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # If more than max_reactants_per_step reactants, it's likely a convergent step
                if len(reactants) > max_reactants_per_step:
                    is_linear = False
                    print(f"Non-linear step detected with {len(reactants)} reactants")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(
        "Synthesis pattern is linear"
        if is_linear
        else "Synthesis pattern is convergent"
    )
    return is_linear
