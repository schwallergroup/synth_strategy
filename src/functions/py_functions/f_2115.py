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
    This function detects if the synthesis follows a linear strategy (as opposed to convergent).
    Linear synthesis is characterized by having mostly 1-2 reactants per step.
    """
    total_reactions = 0
    reactions_with_few_reactants = 0

    def dfs_traverse(node):
        nonlocal total_reactions, reactions_with_few_reactants

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")

            total_reactions += 1
            if len(reactants_smiles) <= 2:
                reactions_with_few_reactants += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # If at least 80% of reactions have few reactants, consider it a linear synthesis
    if total_reactions > 0 and (reactions_with_few_reactants / total_reactions) >= 0.8:
        print(
            f"Linear synthesis detected: {reactions_with_few_reactants}/{total_reactions} reactions have ≤2 reactants"
        )
        return True
    return False
