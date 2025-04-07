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
    Detects if the synthesis follows a linear strategy without convergent steps.
    """
    # Initialize tracking variables
    is_linear = True
    max_reactants_per_step = 0

    def dfs_traverse(node, depth):
        nonlocal is_linear, max_reactants_per_step

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count number of reactants
            num_reactants = len(reactants)
            max_reactants_per_step = max(max_reactants_per_step, num_reactants)

            # If more than 2 complex reactants, it might be convergent
            complex_reactants = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.GetNumAtoms() > 6:  # Arbitrary threshold for "complex"
                    complex_reactants += 1

            if complex_reactants > 2:
                print(
                    f"Potentially convergent step detected at depth {depth} with {complex_reactants} complex reactants"
                )
                is_linear = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route, 0)

    print(f"Linear synthesis strategy detected: {is_linear}")
    return is_linear
