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
    Detects if the synthesis follows a linear strategy without convergent steps.
    """
    # Track the maximum branching factor
    max_branching = 0

    def dfs_traverse(node):
        nonlocal max_branching

        if node["type"] == "reaction":
            if "children" in node:
                # Count number of non-trivial reactants (complex molecules)
                complex_reactants = 0
                for child in node["children"]:
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        # This is a complex intermediate, not a simple starting material
                        complex_reactants += 1

                max_branching = max(max_branching, complex_reactants)

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If max_branching is 1, it's a linear synthesis
    is_linear = max_branching <= 1
    print(f"Maximum branching factor: {max_branching}, Linear synthesis: {is_linear}")
    return is_linear
