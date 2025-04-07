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
    This function detects a synthesis strategy where methoxy groups are preserved
    throughout the synthesis while other functional groups are modified.
    """
    # Track reactions and methoxy presence
    reactions_count = 0
    reactions_with_methoxy = 0

    def dfs_traverse(node):
        nonlocal reactions_count, reactions_with_methoxy

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Count this reaction
            reactions_count += 1

            # Check for methoxy group in product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                methoxy_pattern = Chem.MolFromSmarts("[c][O][C]")
                if product_mol.HasSubstructMatch(methoxy_pattern):
                    reactions_with_methoxy += 1
                    print(f"Found methoxy group in reaction at depth {reactions_count}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if all reactions preserve methoxy groups and we have at least 3 reactions
    result = reactions_count >= 3 and reactions_with_methoxy == reactions_count
    print(
        f"Methoxy preserving synthesis strategy detected: {result} ({reactions_with_methoxy}/{reactions_count} reactions)"
    )
    return result
