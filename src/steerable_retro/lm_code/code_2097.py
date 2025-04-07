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
    Detects if the synthesis follows a linear strategy (each reaction has only one
    non-starting material reactant).
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count non-starting material reactants
            non_starting_material_count = 0

            for reactant in node.get("children", []):
                if reactant["type"] == "mol" and not reactant.get("in_stock", False):
                    non_starting_material_count += 1

            # If more than one non-starting material, it's not linear
            if non_starting_material_count > 1:
                is_linear = False
                print(f"Non-linear step found: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Ensure we have at least 3 reactions to call it a strategy
    if reaction_count < 3:
        is_linear = False

    if is_linear:
        print(f"Linear synthesis strategy detected with {reaction_count} reactions")

    return is_linear
