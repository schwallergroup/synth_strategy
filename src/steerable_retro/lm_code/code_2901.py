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
    This function detects a convergent synthesis strategy with late-stage coupling,
    where multiple fragments are combined in the final steps of the synthesis.
    """
    fragment_count_at_depth = {}
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")

            # Count number of fragments at this depth
            fragment_count = len(reactants)
            if depth in fragment_count_at_depth:
                fragment_count_at_depth[depth] = max(fragment_count_at_depth[depth], fragment_count)
            else:
                fragment_count_at_depth[depth] = fragment_count

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if there are multiple fragments at low depths (late-stage)
    is_convergent = False
    late_stage_threshold = max_depth // 2  # Define late stage as first half of synthesis

    for depth, count in fragment_count_at_depth.items():
        if depth <= late_stage_threshold and count >= 2:
            is_convergent = True
            print(f"Detected convergent synthesis with {count} fragments at depth {depth}")

    return is_convergent
