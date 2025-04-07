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
    This function detects if the synthesis follows a linear strategy (each step has only one non-commercial product as reactant).
    """
    is_linear_synthesis = True

    def dfs_traverse(node):
        nonlocal is_linear_synthesis
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count non-commercial reactants
                non_commercial_count = 0
                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        non_commercial_count += 1

                if non_commercial_count > 1:
                    is_linear_synthesis = False
                    print(
                        f"Non-linear step detected in reaction: {rsmi} with {non_commercial_count} non-commercial reactants"
                    )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Linear synthesis strategy detected: {is_linear_synthesis}")
    return is_linear_synthesis
