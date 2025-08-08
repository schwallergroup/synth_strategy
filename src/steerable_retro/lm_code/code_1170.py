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
    Detects if the synthesis follows a predominantly linear strategy.
    """
    # Count branching points and total reactions in the synthesis route
    branching_points = 0
    total_reactions = 0

    def analyze_structure(node):
        nonlocal branching_points, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1
            # If a reaction node has more than 2 children, it's a branching point
            if len(node.get("children", [])) > 2:
                branching_points += 1

        # Continue DFS traversal
        for child in node.get("children", []):
            analyze_structure(child)

    # Start DFS traversal
    analyze_structure(route)

    # A linear synthesis should have minimal branching and sufficient steps
    is_linear = branching_points <= 1 and total_reactions >= 3
    print(f"Synthesis has {branching_points} branching points and {total_reactions} reactions")
    if is_linear:
        print("Linear synthesis strategy detected")

    return is_linear
