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
    Detects a linear synthesis strategy where each step builds directly on the previous product
    without convergent steps.
    """
    # Track the number of reactants at each step
    reactant_counts = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Count non-empty reactants
                count = sum(1 for r in reactants_smiles if r.strip())

                # Ensure we have enough slots in our list
                while len(reactant_counts) <= depth:
                    reactant_counts.append(0)

                reactant_counts[depth] = count
                print(f"Depth {depth}: {count} reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all steps have at most 2 reactants (linear synthesis)
    # and at least one step has exactly 2 reactants (not just functional group manipulations)
    is_linear = all(count <= 2 for count in reactant_counts) and any(
        count == 2 for count in reactant_counts
    )

    if is_linear:
        print("Detected linear synthesis strategy")
    else:
        print("Not a linear synthesis strategy")

    return is_linear
