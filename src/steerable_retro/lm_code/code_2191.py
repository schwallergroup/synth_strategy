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
    Detects if the synthesis follows a convergent strategy with multiple fragments
    being prepared separately and then combined.
    """
    fragment_paths = []
    max_depth = 0

    def dfs_traverse(node, path=None, depth=0):
        nonlocal fragment_paths, max_depth

        if path is None:
            path = []

        if depth > max_depth:
            max_depth = depth

        # Create a copy of the current path
        current_path = path.copy()

        # Add current node to path if it's a reaction
        if node.get("type") == "reaction":
            # Store the reaction with a unique identifier if available
            if "metadata" in node and "reaction_hash" in node["metadata"]:
                reaction_info = {
                    "node": node,
                    "id": node["metadata"]["reaction_hash"],
                    "depth": depth,
                }
            else:
                # Create a simple identifier based on position in the tree
                reaction_info = {
                    "node": node,
                    "id": f"reaction_at_depth_{depth}_{len(current_path)}",
                    "depth": depth,
                }
            current_path.append(reaction_info)

        # If this is a leaf node (starting material)
        if node.get("type") == "mol" and node.get("in_stock", False):
            fragment_paths.append(current_path)
            return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_path, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have multiple fragment paths
    if len(fragment_paths) < 2:
        print("Not a convergent synthesis: only one fragment path found")
        return False

    print(f"Found {len(fragment_paths)} fragment paths")

    # Check if fragments are combined in the second half of synthesis
    # Find the reaction where fragments are combined
    fragment_combination_depths = []

    for i in range(len(fragment_paths)):
        for j in range(i + 1, len(fragment_paths)):
            path1 = fragment_paths[i]
            path2 = fragment_paths[j]

            # Find common reactions between paths by comparing reaction IDs
            path1_ids = [r["id"] for r in path1]
            path2_ids = [r["id"] for r in path2]
            common_ids = set(path1_ids).intersection(set(path2_ids))

            if common_ids:
                print(f"Found common reactions between path {i} and path {j}")

                # Get the deepest common reaction (first combination point)
                deepest_common = None
                max_common_depth = -1

                for reaction_id in common_ids:
                    # Find the reaction in path1
                    for r in path1:
                        if r["id"] == reaction_id:
                            depth = path1.index(r)
                            if depth > max_common_depth:
                                max_common_depth = depth
                                deepest_common = r

                if deepest_common:
                    fragment_combination_depths.append(max_common_depth)
                    print(f"Deepest common reaction at depth {max_common_depth}")

    if not fragment_combination_depths:
        print("No fragment combination points found")
        return False

    # Calculate the relative position of fragment combination
    avg_combination_depth = sum(fragment_combination_depths) / len(fragment_combination_depths)
    relative_position = avg_combination_depth / max_depth

    # If fragments are combined in the second half of synthesis (relative_position < 0.5)
    # This means the combination happens closer to the target molecule
    is_late_stage = relative_position < 0.5

    print(f"Convergent synthesis detected with {len(fragment_paths)} fragment paths")
    print(f"Fragment combination occurs at relative position: {relative_position:.2f}")
    print(f"Late-stage combination: {is_late_stage}")

    # Always return True if we've detected a convergent synthesis with multiple paths
    # and at least one combination point
    return True
