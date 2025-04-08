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
    This function detects if the synthetic route follows a linear strategy rather than convergent.
    """
    max_branching = 0

    def count_reactants(rsmi):
        reactants = rsmi.split(">")[0].split(".")
        return len(reactants)

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            num_reactants = count_reactants(rsmi)
            max_branching = max(max_branching, num_reactants)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # If max_branching <= 2, it's likely a linear synthesis
    is_linear = max_branching <= 2
    print(f"Maximum branching: {max_branching}, Is linear: {is_linear}")
    return is_linear
