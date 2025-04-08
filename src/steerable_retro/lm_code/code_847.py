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
    rather than a convergent one.
    """
    max_reactants_per_step = 0

    def dfs_traverse(node):
        nonlocal max_reactants_per_step

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count number of reactants in this step
            num_reactants = len(reactants)
            max_reactants_per_step = max(max_reactants_per_step, num_reactants)

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If no step has more than 2 reactants, it's likely a linear synthesis
    is_linear = max_reactants_per_step <= 2
    print(f"Maximum reactants per step: {max_reactants_per_step}, Is linear: {is_linear}")
    return is_linear
