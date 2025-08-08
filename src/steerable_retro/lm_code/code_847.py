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
    This function detects a linear fragment assembly strategy as opposed to convergent.
    """
    # Track reaction depths and branching
    reaction_depths = []
    max_children_per_node = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_children_per_node

        if node["type"] == "reaction":
            reaction_depths.append(depth)

            # Count number of reactants
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                num_reactants = len(reactants)
                max_children_per_node = max(max_children_per_node, num_reactants)

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Analyze reaction pattern
    # Linear synthesis typically has:
    # 1. Monotonically increasing depths
    # 2. Few reactants per step (usually 2)

    is_monotonic = all(
        reaction_depths[i] <= reaction_depths[i + 1] for i in range(len(reaction_depths) - 1)
    )
    few_reactants = max_children_per_node <= 2

    print(f"Reaction depths: {reaction_depths}")
    print(f"Max reactants per step: {max_children_per_node}")
    print(f"Is monotonic: {is_monotonic}")

    return is_monotonic and few_reactants and len(reaction_depths) >= 3
