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
    This function detects if the synthesis follows a linear strategy (no convergent steps).
    """
    is_linear = True

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If there are more than 2 reactants, it might be convergent
            # (allowing for 2 because many reactions have reagents listed as reactants)
            if len(reactants_smiles) > 2:
                is_linear = False
                print("Found potential convergent step with more than 2 reactants")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return is_linear
