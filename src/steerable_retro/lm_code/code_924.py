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
    Detects if the synthesis follows a primarily linear build-up strategy
    rather than a convergent approach.
    """
    # Track reaction depths, branching, and reaction count
    reaction_depths = {}
    max_depth = 0
    total_reactions = 0
    high_branching_reactions = 0

    def dfs_traverse(node, current_depth=0):
        nonlocal max_depth, total_reactions, high_branching_reactions

        if node["type"] == "reaction":
            # Track this reaction's depth
            reaction_depths[current_depth] = reaction_depths.get(current_depth, 0) + 1
            max_depth = max(max_depth, current_depth)
            total_reactions += 1

            # Count reactants (children of reaction node)
            reactant_count = len(node.get("children", []))
            if reactant_count > 2:  # More than 2 reactants suggests convergent synthesis
                high_branching_reactions += 1

            # Process children (reactants) at the same depth level
            for child in node.get("children", []):
                dfs_traverse(child, current_depth + 1)
        else:  # Molecule node
            # Only traverse non-stock molecules (continue synthesis path)
            if not node.get("in_stock", False):
                for child in node.get("children", []):
                    dfs_traverse(child, current_depth)

    # Start traversal
    dfs_traverse(route)

    # Handle empty routes or routes with no reactions
    if total_reactions == 0:
        print("No reactions found in route")
        return False

    # Calculate linearity score
    # In a perfectly linear synthesis, we would have exactly one reaction at each depth
    depth_linearity = sum(1 for depth, count in reaction_depths.items() if count == 1) / max(
        1, len(reaction_depths)
    )

    # Calculate branching score (lower is more linear)
    branching_score = high_branching_reactions / max(1, total_reactions)

    # Combined linearity metric
    linearity_score = depth_linearity * (1 - branching_score)

    print(f"Linearity score: {linearity_score:.2f}, reaction depths: {reaction_depths}")
    print(f"Depth linearity: {depth_linearity:.2f}, Branching score: {branching_score:.2f}")

    # If most depths have only one reaction and branching is low, it's a linear strategy
    return linearity_score > 0.6
