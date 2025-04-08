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
    (each step builds on a single precursor, not convergent)
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # If a reaction has more than one non-reagent reactant, it's convergent
            # We'll consider small molecules (less than 10 atoms) as potential reagents
            significant_reactants = 0
            for reactant in reactants:
                mol = Chem.MolFromSmiles(reactant)
                if (
                    mol and mol.GetNumAtoms() >= 10
                ):  # Arbitrary threshold for "significant" molecule
                    significant_reactants += 1

            if significant_reactants > 1:
                is_linear = False
                print(f"Found convergent step with {significant_reactants} significant reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Linear synthesis strategy detected: {is_linear}")
    return is_linear
