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
    This function detects if the synthetic route follows a linear strategy
    (each reaction typically has 2 reactants, one being the product of the previous step).
    """
    # For a linear synthesis, we expect most reactions to have 2 reactants
    linear_strategy = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal linear_strategy, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                reactant_list = reactants_smiles.split(".")

                # If more than 2 reactants, it might not be a linear synthesis
                if len(reactant_list) > 2:
                    linear_strategy = False
                    print(
                        f"Non-linear step detected with {len(reactant_list)} reactants"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Need at least 3 reactions to confirm linear strategy
    if reaction_count < 3:
        return False

    return linear_strategy
