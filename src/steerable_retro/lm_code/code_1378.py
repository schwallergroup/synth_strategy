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
    This function detects if the synthesis uses a Grignard reaction for C-C bond formation.
    """
    uses_grignard = False

    def dfs_traverse(node):
        nonlocal uses_grignard

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]

            # Check if any reactant contains Mg (characteristic of Grignard reagents)
            if "Mg" in reactants_part:
                print("Grignard reagent detected in reactants")
                uses_grignard = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return uses_grignard
