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
    This function detects if the synthesis follows a linear pattern where
    each step builds upon a single product from the previous step.
    """
    # Track the number of reactants at each step
    reaction_counts = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]

                # Count the number of reactants (separated by '.')
                reactant_count = len(reactants_part.split("."))
                reaction_counts.append((depth, reactant_count))
                print(f"Reaction at depth {depth} has {reactant_count} reactants")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # A linear synthesis typically has 2 reactants per step (the main product from previous step + one reagent)
    # We'll consider it linear if most reactions have 2 or fewer reactants
    linear_reactions = [count for depth, count in reaction_counts if count <= 2]
    is_linear = (
        len(linear_reactions) >= len(reaction_counts) * 0.75
    )  # At least 75% of reactions are linear

    print(f"Linear synthesis strategy detected: {is_linear}")
    return is_linear
