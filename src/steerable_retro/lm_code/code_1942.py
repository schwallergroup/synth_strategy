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
    This function detects if the synthetic route follows a linear synthesis strategy
    (each step builds upon a single precursor molecule).
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If more than 2 reactants, it's not a linear synthesis
            if len(reactants_smiles) > 2:
                is_linear = False
                print(f"More than 2 reactants found: {len(reactants_smiles)}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return is_linear
