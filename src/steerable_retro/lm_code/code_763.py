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
    Detects if the synthesis has a branched structure with multiple starting materials.
    """
    # Count the number of leaf nodes (starting materials)
    leaf_count = 0

    def dfs_traverse(node):
        nonlocal leaf_count

        if node["type"] == "mol" and node.get("in_stock", False):
            leaf_count += 1
            return

        if not node.get("children", []):
            # This is a leaf node but not marked as in_stock
            leaf_count += 1
            return

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # A branched synthesis typically has multiple starting materials
    is_branched = leaf_count >= 3
    if is_branched:
        print(f"Detected branched synthesis with {leaf_count} starting materials")

    return is_branched
