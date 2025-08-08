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
    This function detects if the route follows a linear synthesis strategy
    (as opposed to convergent synthesis with multiple fragments).

    A synthesis is considered linear if:
    1. Most reactions (>80%) have 2 or fewer reactants
    2. No reaction has more than 3 reactants
    """
    max_fragments_per_reaction = 0
    total_reactions = 0
    reactions_with_many_fragments = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_fragments_per_reaction, total_reactions, reactions_with_many_fragments

        if node["type"] == "reaction":
            try:
                if "rsmi" in node["metadata"]:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    num_fragments = len(reactants)
                    max_fragments_per_reaction = max(max_fragments_per_reaction, num_fragments)
                    total_reactions += 1

                    if num_fragments > 2:
                        reactions_with_many_fragments += 1

                    print(f"Reaction at depth {depth} has {num_fragments} fragments")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # No reactions found
    if total_reactions == 0:
        return True

    # Calculate percentage of reactions with more than 2 fragments
    percentage_many_fragments = (reactions_with_many_fragments / total_reactions) * 100
    print(f"Max fragments per reaction: {max_fragments_per_reaction}")
    print(f"Percentage of reactions with >2 fragments: {percentage_many_fragments:.1f}%")

    # A synthesis is linear if:
    # 1. No reaction has more than 3 fragments
    # 2. Less than or equal to 20% of reactions have more than 2 fragments
    return max_fragments_per_reaction <= 3 and percentage_many_fragments <= 20
