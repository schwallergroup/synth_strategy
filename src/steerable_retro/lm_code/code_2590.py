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
    Detects synthesis routes that use a linear fragment coupling strategy
    (as opposed to convergent synthesis).
    """
    # Track the number of fragments combined at each step
    fragment_counts = []
    # Track the depth of each reaction for structure analysis
    reaction_depths = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal fragment_counts, reaction_depths, max_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")

            # Count the number of fragments (reactants)
            num_fragments = len([r for r in reactants if r])

            # Only consider coupling reactions (2 or more fragments)
            if num_fragments >= 2:
                fragment_counts.append(num_fragments)
                reaction_depths[len(fragment_counts) - 1] = depth
                print(f"Found coupling reaction with {num_fragments} fragments at depth {depth}")
            else:
                print(
                    f"Skipping non-coupling reaction with {num_fragments} fragment at depth {depth}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have enough coupling reactions
    if len(fragment_counts) >= 2:
        # Check if most coupling reactions involve only 2 fragments (linear synthesis)
        two_fragment_reactions = sum(1 for count in fragment_counts if count == 2)
        linear_ratio = two_fragment_reactions / len(fragment_counts)
        print(f"Linear fragment coupling ratio: {linear_ratio}")

        # Check if the reaction structure is linear
        # In a linear structure, reactions should occur at different depths
        unique_depths = len(set(reaction_depths.values()))
        depth_ratio = unique_depths / len(reaction_depths) if reaction_depths else 0
        print(f"Depth distribution ratio: {depth_ratio}")

        # A route is linear if:
        # 1. At least 60% of coupling reactions are binary (2 fragments)
        # 2. Reactions occur at different depths (indicating sequential rather than parallel)
        is_linear = linear_ratio >= 0.6 and depth_ratio >= 0.7

        return is_linear

    print("Not enough coupling reactions to determine strategy")
    return False
