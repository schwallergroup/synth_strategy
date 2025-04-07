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
    This function detects if the synthetic route follows a linear synthesis strategy
    without convergent steps.
    """
    is_linear = True
    reaction_counts = {}

    def dfs_traverse(node, depth=0):
        nonlocal is_linear, reaction_counts

        if node["type"] == "reaction":
            # Track reactions at each depth level
            reaction_counts[depth] = reaction_counts.get(depth, 0) + 1

            # Check if this reaction has multiple reactants that are products of previous reactions
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")

                # Count non-starting material reactants
                non_starting_material_reactants = 0
                for child in node.get("children", []):
                    if child["type"] == "mol" and not child.get("in_stock", False):
                        non_starting_material_reactants += 1

                if non_starting_material_reactants > 1:
                    print(f"Found potential convergent step at depth {depth}: {rsmi}")
                    is_linear = False

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have more than one reaction at any depth level
    # This indicates branching in the synthesis
    for depth, count in reaction_counts.items():
        if count > 1:
            print(
                f"Found multiple reactions at depth {depth}, indicating non-linear synthesis"
            )
            is_linear = False

    return is_linear
