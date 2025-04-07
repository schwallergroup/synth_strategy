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
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # In a linear synthesis, we typically have 1-2 reactants per step
            # More than 2 reactants often indicates a convergent approach
            if len(reactants_smiles) > 2:
                is_linear = False
                print(f"Non-linear step detected with {len(reactants_smiles)} reactants")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    print(f"Has linear synthesis strategy: {is_linear}")
    return is_linear
