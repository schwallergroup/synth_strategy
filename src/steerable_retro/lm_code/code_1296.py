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
    Detects if the synthesis follows a linear build-up strategy rather than
    a convergent approach by analyzing the reaction tree structure.
    """
    # Track branching factor at each depth
    depth_to_children = {}
    max_branching = 0

    def dfs_traverse(node, current_depth=0):
        nonlocal max_branching

        if node["type"] == "reaction":
            depth = node["metadata"].get("depth", current_depth)
            children = node.get("children", [])
            child_count = len(children)

            if depth not in depth_to_children:
                depth_to_children[depth] = child_count
            else:
                depth_to_children[depth] += child_count

            max_branching = max(max_branching, child_count)

            for child in children:
                dfs_traverse(child, depth + 1)
        else:
            for child in node.get("children", []):
                dfs_traverse(child, current_depth)

    dfs_traverse(route)

    # If max branching factor is <= 2, it's likely a linear synthesis
    # (allowing for 2 in case of protection groups or similar)
    is_linear = max_branching <= 2

    print(
        f"Synthesis has max branching factor of {max_branching}, classified as {'linear' if is_linear else 'convergent'}"
    )
    return is_linear
