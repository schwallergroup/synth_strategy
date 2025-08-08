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
    This function detects if the synthesis follows a linear strategy with no branching.
    """
    is_linear = True
    max_children_per_node = 0

    def dfs_traverse(node):
        nonlocal is_linear, max_children_per_node

        # Count children
        num_children = len(node.get("children", []))
        max_children_per_node = max(max_children_per_node, num_children)

        # If a node has more than 2 children, it's not linear
        if num_children > 2:
            is_linear = False
            print("Found node with more than 2 children, not a linear synthesis")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # A linear synthesis should have at most 2 children per node
    # (1 for the next reaction, 1 for a potential reagent)
    return is_linear and max_children_per_node <= 2
