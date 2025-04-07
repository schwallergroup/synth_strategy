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
    This function detects if the synthesis follows a linear path without convergent steps.
    A linear synthesis has at least 3 reactions and each reaction has at most one non-starting material reactant.
    """
    is_linear = True
    reaction_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, reaction_count

        if node["type"] == "reaction":
            reaction_count += 1

            # Check if this reaction has more than one non-starting material reactant
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")

                # Find the children nodes that correspond to reactants
                children = node.get("children", [])

                # Count non-starting material reactants
                complex_reactants = 0
                for child in children:
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        complex_reactants += 1

                if complex_reactants > 1:
                    is_linear = False
                    print(f"Convergent step detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # A synthesis with at least 3 reactions that is linear
    return is_linear and reaction_count >= 3
