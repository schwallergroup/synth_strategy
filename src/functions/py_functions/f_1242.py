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
    Detects if the synthesis follows a linear strategy rather than a convergent one.
    Linear synthesis typically has one main reactant and one smaller fragment in each step.
    """
    is_linear = True
    step_count = 0

    def dfs_traverse(node):
        nonlocal is_linear, step_count

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            step_count += 1

            # If any step has more than 2 reactants, it might be convergent
            if len(reactants) > 2:
                is_linear = False
                print(
                    f"Detected potential convergent step with {len(reactants)} reactants at depth {node['metadata'].get('depth', 'unknown')}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Must have at least 3 steps to be considered a meaningful linear synthesis
    return is_linear and step_count >= 3
