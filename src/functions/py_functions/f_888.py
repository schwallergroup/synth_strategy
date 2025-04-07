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
    Detects if the synthesis follows a linear pattern where each step
    builds upon a single product from the previous step.

    In a linear synthesis:
    1. Each reaction should have exactly one product molecule
    2. Each intermediate molecule (not a starting material) should be used in exactly one reaction
    3. The synthesis forms a single chain without branches
    """
    is_linear = True

    # Track visited nodes to avoid cycles
    visited = set()

    def dfs_traverse(node, depth=0):
        nonlocal is_linear

        # Skip if already visited
        node_id = id(node)
        if node_id in visited:
            return
        visited.add(node_id)

        if node["type"] == "reaction":
            # For a reaction node, check if it produces exactly one product molecule
            # In a retrosynthetic tree, we need to check the reaction SMILES
            # to determine the number of products in the forward direction
            try:
                rsmi = node["metadata"]["rsmi"]
                products = rsmi.split(">")[-1].split(".")
                product_count = len(products)

                if product_count != 1:
                    is_linear = False
                    print(
                        f"Non-linear pattern detected at depth {depth}: reaction has {product_count} product molecules in forward direction"
                    )
                    print(f"Reaction SMILES: {rsmi}")
            except (KeyError, IndexError) as e:
                # If we can't get the reaction SMILES, assume it's linear and continue
                print(
                    f"Warning at depth {depth}: Could not analyze reaction SMILES - {str(e)}"
                )

        elif node["type"] == "mol":
            # Skip checking for the target molecule (depth 0) and starting materials
            if depth > 0 and not node.get("in_stock", False):
                # Count how many reactions this molecule participates in
                reaction_children = [
                    child
                    for child in node.get("children", [])
                    if child["type"] == "reaction"
                ]

                # In a linear synthesis, each intermediate molecule should be used in exactly one reaction
                if len(reaction_children) != 1:
                    is_linear = False
                    print(
                        f"Non-linear pattern detected at depth {depth}: molecule feeds into {len(reaction_children)} reactions"
                    )
                    print(f"Molecule SMILES: {node.get('smiles', 'Not available')}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return is_linear
