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
    This function detects if the synthesis follows a linear strategy without convergent steps.
    """
    is_linear = True
    step_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, step_count

        if node["type"] == "reaction":
            step_count += 1
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # If more than 2 reactants, it might be convergent
            if len(reactants_smiles) > 2:
                is_linear = False
                print(
                    f"Convergent step detected with {len(reactants_smiles)} reactants"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have a linear synthesis with at least 3 steps
    is_linear_synthesis = is_linear and step_count >= 3

    print(
        f"Linear synthesis strategy detected: {is_linear_synthesis} (steps: {step_count})"
    )
    return is_linear_synthesis
