#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects if the synthesis follows a linear strategy where each step builds
    on the previous with no convergent steps.

    A linear synthesis has a tree structure with no branching except possibly
    at the final step (depth 0).
    """
    # Track the tree structure to detect branching
    mol_nodes_by_depth = {}

    def dfs_traverse(node, depth=0, current_path=None):
        if current_path is None:
            current_path = []

        if node["type"] == "mol":
            # Skip in-stock molecules (starting materials)
            if node.get("in_stock", False):
                return

            # Store molecule at its depth
            if depth not in mol_nodes_by_depth:
                mol_nodes_by_depth[depth] = []
            mol_nodes_by_depth[depth].append(node["smiles"])

        # Add current node to path and traverse children
        new_path = current_path + [node]
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, new_path)

    # Start traversal from the root
    dfs_traverse(route)

    # If no non-stock molecules were found, return False
    if not mol_nodes_by_depth:
        print("Linear synthesis strategy detected: False (no non-stock molecules found)")
        return False

    # Sort depths (lower depth = later in synthesis, closer to final product)
    sorted_depths = sorted(mol_nodes_by_depth.keys())

    # A linear synthesis should have exactly one molecule at each depth
    # except possibly at the lowest depth (final product)
    linear = True
    if sorted_depths:
        min_depth = min(sorted_depths)
        for depth in sorted_depths:
            if depth != min_depth and len(mol_nodes_by_depth[depth]) > 1:
                linear = False
                print(
                    f"Non-linear branching detected at depth {depth} with {len(mol_nodes_by_depth[depth])} molecules"
                )
                break

    print(f"Linear synthesis strategy detected: {linear}")
    return linear
