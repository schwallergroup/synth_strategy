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
    This function detects if the synthesis follows a linear build-up strategy
    where each step builds on the previous product without convergent steps.
    """
    is_linear = True

    def dfs_traverse(node):
        nonlocal is_linear

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count non-starting material reactants
                non_starting_material_count = 0

                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        non_starting_material_count += 1

                # If more than one non-starting material reactant, it's convergent
                if non_starting_material_count > 1:
                    is_linear = False
                    print("Convergent step detected - not a linear build-up strategy")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if is_linear:
        print("Linear build-up strategy confirmed")

    return is_linear
