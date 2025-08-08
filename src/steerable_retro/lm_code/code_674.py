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
    Detects a linear synthesis pathway with sequential transformations.
    """
    # Track reaction depths and branching
    reaction_depths = set()
    branching_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal reaction_depths, branching_detected

        if node["type"] == "reaction":
            reaction_depths.add(depth)

            # Check for branching (multiple children for a reaction node)
            reaction_children = [c for c in node.get("children", []) if c["type"] == "reaction"]
            if len(reaction_children) > 1:
                branching_detected = True
                print(f"Branching detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Linear synthesis has sequential reactions without branching
    strategy_present = len(reaction_depths) >= 3 and not branching_detected

    if strategy_present:
        print(f"Linear synthesis detected with {len(reaction_depths)} sequential reactions")
    else:
        print("Linear synthesis strategy not detected")

    return strategy_present
