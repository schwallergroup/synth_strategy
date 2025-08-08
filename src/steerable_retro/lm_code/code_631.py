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
    This function detects a strategy where the synthesis starts with a convergent step
    (multiple reactants) followed by linear build-up (single reactant per step).
    """
    # Store paths from target to starting materials
    all_paths = []

    def dfs_traverse(node, depth=0, path=None):
        if path is None:
            path = []

        current_path = path.copy()

        # If this is a reaction node, classify it
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]

            # Count number of reactants
            reactant_count = len([r for r in reactants_part.split(".") if r])

            reaction_type = "convergent" if reactant_count > 1 else "linear"
            current_path.append((depth, reaction_type))

            print(f"Found {reaction_type} step with {reactant_count} reactants at depth {depth}")

        # If this is a leaf node (starting material), save the path
        if not node.get("children", []):
            if current_path:  # Only save if the path contains reactions
                all_paths.append(current_path)
                print(f"Saved path with {len(current_path)} reactions: {current_path}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal from the target molecule
    dfs_traverse(route)

    # Check if any path has the pattern we're looking for
    for path in all_paths:
        if not path:  # Skip empty paths
            continue

        # In retrosynthesis, we traverse from target to starting materials
        # So the first reactions in the path are late-stage, and last ones are early-stage

        # Find the earliest convergent step (highest depth)
        convergent_indices = [
            i for i, (_, reaction_type) in enumerate(path) if reaction_type == "convergent"
        ]

        if not convergent_indices:
            continue  # No convergent steps in this path

        earliest_convergent_idx = max(convergent_indices, key=lambda i: path[i][0])
        earliest_convergent_depth = path[earliest_convergent_idx][0]

        # Check if all steps after the earliest convergent step (toward the target) are linear
        later_steps = [
            reaction_type for depth, reaction_type in path if depth < earliest_convergent_depth
        ]

        if all(step == "linear" for step in later_steps):
            print(
                f"Detected convergent step at depth {earliest_convergent_depth} followed by linear build-up"
            )
            return True

    print("Linear build-up after convergent step strategy not detected")
    return False
