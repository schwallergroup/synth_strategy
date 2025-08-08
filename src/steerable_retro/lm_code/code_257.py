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
    This function detects a linear synthesis strategy that builds complexity
    through sequential transformations rather than convergent assembly.
    """
    # Track reaction counts and path structure
    reaction_count = 0
    branch_points = 0
    max_path_length = 0
    paths = []

    # Track the structure of the synthesis tree
    def dfs_traverse(node, current_depth=0, current_path=None):
        nonlocal reaction_count, branch_points, max_path_length

        if current_path is None:
            current_path = []

        # Update max path length
        max_path_length = max(max_path_length, current_depth)

        # Add current node to path
        current_path.append(node)

        if node["type"] == "reaction":
            reaction_count += 1

            # Count children that are molecules and not in stock
            non_stock_children = [
                child
                for child in node.get("children", [])
                if child["type"] == "mol" and not child.get("in_stock", False)
            ]

            # If more than one non-stock child, this is a potential branch point
            # But we need to check if they're part of the same linear sequence
            if len(non_stock_children) > 1:
                # Check if this is a true convergent step or just a reaction with multiple products
                # For true convergent synthesis, the non-stock children should lead to different reaction paths
                child_reaction_counts = []
                for child in non_stock_children:
                    child_reactions = 0

                    def count_child_reactions(n):
                        nonlocal child_reactions
                        if n["type"] == "reaction":
                            child_reactions += 1
                        for c in n.get("children", []):
                            count_child_reactions(c)

                    count_child_reactions(child)
                    child_reaction_counts.append(child_reactions)

                # If multiple children have significant reaction paths, this is a true branch point
                significant_paths = sum(1 for count in child_reaction_counts if count > 0)
                if significant_paths > 1:
                    branch_points += 1

        # If this is a leaf node (no children), save the path
        if not node.get("children", []):
            paths.append(list(current_path))
        else:
            # Traverse children
            for child in node.get("children", []):
                dfs_traverse(child, current_depth + 1, list(current_path))

    # Start traversal
    dfs_traverse(route)

    # Find the longest path with reactions
    longest_path_reactions = 0
    for path in paths:
        path_reactions = sum(1 for node in path if node["type"] == "reaction")
        longest_path_reactions = max(longest_path_reactions, path_reactions)

    # Calculate linearity ratio - lower means more linear
    linearity_ratio = (branch_points + 1) / (reaction_count + 1) if reaction_count > 0 else 1

    # Calculate path dominance - higher means more linear
    path_dominance = longest_path_reactions / reaction_count if reaction_count > 0 else 0

    # Linear synthesis should have:
    # 1. A reasonable ratio of branch points to reactions
    # 2. At least 2 reaction steps
    # 3. A dominant linear path containing most reactions
    strategy_present = linearity_ratio <= 0.7 and reaction_count >= 2 and path_dominance >= 0.5

    print(f"Strategy detection results:")
    print(f"  Reaction count: {reaction_count}")
    print(f"  Branch points: {branch_points}")
    print(f"  Max path length: {max_path_length}")
    print(f"  Linearity ratio: {linearity_ratio:.2f}")
    print(f"  Longest path reactions: {longest_path_reactions}")
    print(f"  Path dominance: {path_dominance:.2f}")
    print(f"  Linear complexity buildup: {strategy_present}")

    return strategy_present
