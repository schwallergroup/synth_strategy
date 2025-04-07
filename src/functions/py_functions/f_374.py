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
    This function detects if the synthesis follows a linear strategy where each step
    builds on a single previous intermediate (as opposed to convergent synthesis).
    """
    is_linear = True

    def count_non_reagent_reactants(rsmi):
        """Count reactants that are likely not just reagents based on complexity"""
        reactants = rsmi.split(">")[0].split(".")
        # Filter out small molecules that are likely reagents
        significant_reactants = [
            r for r in reactants if len(r) > 10 and ("[c]" in r or "[C]" in r)
        ]
        return len(significant_reactants)

    def dfs_traverse(node):
        nonlocal is_linear

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            # If any step has more than one significant reactant, it's not a linear synthesis
            if count_non_reagent_reactants(rsmi) > 1:
                is_linear = False
                print("Found a convergent step with multiple significant reactants")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Linear synthesis detection result: {is_linear}")
    return is_linear
