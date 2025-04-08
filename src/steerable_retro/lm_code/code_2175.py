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
    Detects if the synthesis follows a convergent pattern where multiple fragments
    are combined in a non-linear fashion.
    """
    # Track branching in the synthesis tree
    branch_points = 0
    multi_component_reactions = 0

    def dfs_traverse(node):
        nonlocal branch_points, multi_component_reactions

        # Count children that are reactions
        reaction_children = [
            child for child in node.get("children", []) if child["type"] == "reaction"
        ]

        if len(reaction_children) > 1:
            branch_points += 1
            print(f"Detected branch point with {len(reaction_children)} reaction paths")

        # Check for multi-component reactions (3+ reactants)
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            reactant_count = len(reactants_part.split("."))

            if reactant_count >= 2:
                multi_component_reactions += 1
                print(f"Detected multi-component reaction with {reactant_count} reactants")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    # Consider it convergent if there are branch points or multi-component reactions
    return branch_points > 0 or multi_component_reactions >= 2
