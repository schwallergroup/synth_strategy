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
    adds one fragment to build the molecule (as opposed to convergent synthesis).
    """
    # In a linear synthesis, most reactions have only one product-contributing reactant
    linear_steps = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal linear_steps, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count significant reactants (those with more than 6 atoms)
                significant_reactants = 0
                for r in reactants:
                    if Chem.MolFromSmiles(r) is not None:
                        if Chem.MolFromSmiles(r).GetNumAtoms() > 6:
                            significant_reactants += 1

                if significant_reactants <= 1:
                    linear_steps += 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # If more than 75% of steps are linear, consider it a linear synthesis
    is_linear = (total_steps > 0) and (linear_steps / total_steps >= 0.75)
    if is_linear:
        print(f"Linear synthesis strategy detected: {linear_steps}/{total_steps} steps are linear")

    return is_linear
