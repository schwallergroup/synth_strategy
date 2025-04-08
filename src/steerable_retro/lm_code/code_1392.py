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
    Detects if the synthesis follows a linear strategy (no convergent steps).
    Linear synthesis typically has 1-2 reactants per step.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If more than 2 reactants, likely not linear
            if len(reactants_smiles) > 2:
                is_linear = False
                print(f"Non-linear step detected with {len(reactants_smiles)} reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if is_linear:
        print("Linear synthesis strategy detected")

    return is_linear
