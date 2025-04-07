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
    This function detects a linear fragment assembly strategy (as opposed to convergent).
    """
    reaction_counts = 0
    multi_reactant_counts = 0

    def dfs_traverse(node):
        nonlocal reaction_counts, multi_reactant_counts

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            reaction_counts += 1
            if len(reactants) >= 2:
                multi_reactant_counts += 1
                print(f"Detected multi-reactant step: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis has most steps with single reactant or simple transformations
    # If more than 70% of reactions have multiple reactants, it's likely convergent
    is_linear = reaction_counts > 0 and (multi_reactant_counts / reaction_counts) < 0.7
    print(
        f"Linear fragment assembly detected: {is_linear} (multi-reactant steps: {multi_reactant_counts}/{reaction_counts})"
    )
    return is_linear
