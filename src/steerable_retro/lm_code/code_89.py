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
    This function detects if the synthesis follows a linear strategy (as opposed to convergent).
    Linear synthesis typically has fewer reactants per step.
    """
    reactant_counts = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                reactant_counts.append(len(reactants))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Handle edge case of empty reactant_counts
    if not reactant_counts:
        print("No reaction nodes found in the route.")
        return False

    # Calculate percentage of reactions with 2 or fewer reactants
    percentage = sum(count <= 2 for count in reactant_counts) / len(reactant_counts)

    # If most reactions have 2 or fewer reactants, it's likely a linear synthesis
    # Using 60% threshold instead of 70%
    is_linear = percentage >= 0.6

    print(f"Linear synthesis strategy detected: {is_linear}")
    print(f"Reactant counts per step: {reactant_counts}")
    print(f"Percentage of steps with ≤2 reactants: {percentage:.2f} (threshold: 0.6)")

    return is_linear
