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
    (each reaction typically has 2 reactants).
    """
    reaction_count = 0
    linear_reactions = 0

    def dfs_traverse(node):
        nonlocal reaction_count, linear_reactions

        if node["type"] == "reaction":
            # Extract reactants from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            reaction_count += 1

            # Linear reactions typically have 2 reactants
            if len(reactants) == 2:
                linear_reactions += 1
                print(f"Found linear reaction with 2 reactants")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Consider it a linear synthesis if most reactions have 2 reactants
    if reaction_count > 0:
        return linear_reactions / reaction_count >= 0.75

    return False
