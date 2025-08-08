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
    This function detects a linear synthesis with sequential functionalization
    of an aromatic core.
    """
    # Track reactions by depth
    reactions_by_depth = {}
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if depth not in reactions_by_depth:
                reactions_by_depth[depth] = []

            reactions_by_depth[depth].append(node["metadata"]["rsmi"])

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if we have a linear sequence (one reaction per depth)
    is_linear = all(len(reactions) == 1 for reactions in reactions_by_depth.values())

    # Check if we have at least 3 steps
    has_multiple_steps = len(reactions_by_depth) >= 3

    if is_linear and has_multiple_steps:
        print(f"Linear synthesis with {len(reactions_by_depth)} sequential steps detected")
        return True
    return False
