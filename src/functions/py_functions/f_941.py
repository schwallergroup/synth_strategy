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
    Detects convergent synthesis by checking if multiple fragments are combined
    in a non-linear fashion.
    """
    # Track reaction depths and number of reactants
    reaction_data = []

    def dfs_traverse(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            reaction_data.append({"depth": depth, "num_reactants": len(reactants)})

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check for convergent synthesis patterns
    # 1. Multiple reactions with more than one reactant
    multi_reactant_reactions = [r for r in reaction_data if r["num_reactants"] > 1]

    # 2. Reactions at different depths (non-linear)
    unique_depths = len(set(r["depth"] for r in reaction_data))

    is_convergent = len(multi_reactant_reactions) >= 2 and unique_depths >= 2

    if is_convergent:
        print("Convergent synthesis detected with multiple fragment combinations")

    return is_convergent
