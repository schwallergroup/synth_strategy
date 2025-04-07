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
    Detects if the synthesis involves a late-stage coupling of multiple fragments
    (3 or more) to form the final product.
    """
    found_multi_fragment = False

    def dfs_traverse(node, depth=0):
        nonlocal found_multi_fragment

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_str = rsmi.split(">")[0]

                # Count the number of distinct reactant fragments
                reactants = reactants_str.split(".")
                if len(reactants) >= 3:
                    print(
                        f"Found multi-fragment coupling with {len(reactants)} fragments at depth {depth}"
                    )
                    found_multi_fragment = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return found_multi_fragment
