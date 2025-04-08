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
    # Track if we've found a stereocenter in an early step
    early_stereocenter = False
    linear_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal early_stereocenter, linear_steps

        if node["type"] == "reaction":
            linear_steps += 1

            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product_smiles = rsmi.split(">")[-1]

                # Check for stereochemistry in early steps (depth >= 3)
                if depth >= 3 and "@" in product_smiles:
                    print(f"Detected stereocenter in early step at depth {depth}")
                    early_stereocenter = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # Return True if we have a linear synthesis (4+ steps) with early stereocenter
    return linear_steps >= 4 and early_stereocenter
