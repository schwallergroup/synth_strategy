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


def main(route, min_depth=5):
    """
    This function detects if the synthesis follows a linear strategy with at least a minimum depth.
    """
    max_depth = 0
    branching_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, branching_detected

        max_depth = max(max_depth, depth)

        # Check for branching (more than one child for a molecule node)
        if node["type"] == "mol" and len(node.get("children", [])) > 1:
            branching_detected = True
            print(f"Branching detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = max_depth >= min_depth and not branching_detected
    print(f"Linear synthesis strategy (depth {max_depth} >= {min_depth}, no branching): {result}")
    return result
