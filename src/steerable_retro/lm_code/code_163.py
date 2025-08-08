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
    This function detects if the synthesis follows a linear strategy with at least 3 steps.
    """
    reaction_depths = set()
    max_children_per_node = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_children_per_node

        if node["type"] == "reaction":
            reaction_depths.add(depth)

        children_count = len(node.get("children", []))
        max_children_per_node = max(max_children_per_node, children_count)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Linear synthesis has consecutive depths and no more than 2 children per node
    is_linear = len(reaction_depths) >= 3 and max_children_per_node <= 2

    print(f"Linear synthesis strategy: {is_linear}")
    return is_linear
