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
    This function detects if the synthetic route follows a linear build-up strategy
    rather than a convergent approach.
    """
    # Track the maximum number of reactants in any single step
    max_reactants = 0

    def dfs_traverse(node):
        nonlocal max_reactants

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            # Count non-empty reactants
            num_reactants = sum(1 for r in reactants_smiles if r.strip())
            max_reactants = max(max_reactants, num_reactants)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Linear build-up typically has no more than 2 reactants per step
    # (one main molecule and one reagent)
    print(f"Maximum reactants in any step: {max_reactants}")
    return max_reactants <= 2
