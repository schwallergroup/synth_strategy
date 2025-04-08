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
    Detects if fragment coupling occurs in the second half of the synthesis
    """
    max_depth = 0
    fragment_coupling_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal max_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if this is a fragment coupling (multiple reactants to single product)
            if len([r for r in reactants_smiles if r]) > 1:
                fragment_coupling_depths.append(depth)
                print(f"Detected fragment coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if any fragment coupling occurs in the second half of the synthesis
    if not fragment_coupling_depths:
        return False

    half_depth = max_depth / 2
    late_couplings = [d for d in fragment_coupling_depths if d < half_depth]

    print(f"Fragment couplings: {fragment_coupling_depths}, Max depth: {max_depth}")
    print(f"Late stage couplings: {late_couplings}")

    return len(late_couplings) > 0
