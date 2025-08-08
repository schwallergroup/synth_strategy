#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    This function detects if the route follows a linear synthesis strategy rather than convergent.

    In a linear synthesis:
    1. Most reaction nodes have only one non-starting-material child
    2. Most molecule nodes (except starting materials) have only one parent reaction
    3. The overall structure resembles a "backbone" with minimal branching
    """
    # Track branching statistics
    linear_reactions = 0
    total_reactions = 0
    branching_molecules = 0
    total_non_starting_molecules = 0

    # Track the maximum branching factor in the route
    max_branching = 0

    def dfs_traverse(node, depth=0):
        nonlocal linear_reactions, total_reactions, branching_molecules, total_non_starting_molecules, max_branching

        if node["type"] == "reaction":
            total_reactions += 1

            # Count non-starting material children
            non_starting_children = [
                child
                for child in node.get("children", [])
                if child["type"] == "mol" and not child.get("in_stock", False)
            ]

            # In linear synthesis, a reaction should have exactly one non-starting material child
            if len(non_starting_children) == 1:
                linear_reactions += 1

            # Track maximum branching at any reaction
            current_branching = len(non_starting_children)
            max_branching = max(max_branching, current_branching)

        elif node["type"] == "mol" and not node.get("in_stock", False):
            total_non_starting_molecules += 1

            # Count reaction children for this molecule
            reaction_children = [
                child for child in node.get("children", []) if child["type"] == "reaction"
            ]

            # In a convergent synthesis, molecules often branch into multiple reactions
            if len(reaction_children) > 1:
                branching_molecules += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Calculate linearity metrics
    reaction_linearity = linear_reactions / total_reactions if total_reactions > 0 else 0
    molecule_non_branching = 1 - (
        branching_molecules / total_non_starting_molecules
        if total_non_starting_molecules > 0
        else 0
    )

    print(f"Reaction linearity: {reaction_linearity:.2f} ({linear_reactions}/{total_reactions})")
    print(
        f"Molecule non-branching: {molecule_non_branching:.2f} ({total_non_starting_molecules-branching_molecules}/{total_non_starting_molecules})"
    )
    print(f"Maximum branching factor: {max_branching}")

    # A route is considered linear if:
    # 1. At least 80% of reactions have only one non-starting material child
    # 2. Few molecules branch into multiple reactions
    # 3. The maximum branching factor is low
    is_linear = reaction_linearity >= 0.8 and molecule_non_branching >= 0.8 and max_branching <= 2

    return is_linear
