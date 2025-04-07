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
    Detects if the synthesis follows a linear (non-convergent) approach.
    """
    is_linear = True
    max_reactants_per_step = 2  # Allow up to 2 reactants per step for linear synthesis

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If more than max_reactants_per_step, it's likely convergent
            if len(reactants) > max_reactants_per_step:
                is_linear = False
                print(
                    f"Found {len(reactants)} reactants in a step, synthesis not linear"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
