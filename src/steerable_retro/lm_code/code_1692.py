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
    This function detects if the synthesis follows a linear pattern without convergent steps.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if rsmi:
                reactants = rsmi.split(">")[0].split(".")
                # Count non-empty reactants
                reactant_count = sum(1 for r in reactants if r.strip())

                # If more than 2 reactants, it's likely a convergent step
                if reactant_count > 2:
                    print(f"Convergent step detected with {reactant_count} reactants")
                    is_linear = False

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return is_linear
