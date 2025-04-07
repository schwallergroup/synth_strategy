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
    This function detects if a tert-butyl group is maintained throughout the synthesis.
    """
    tert_butyl_at_depths = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal tert_butyl_at_depths, max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None:
                tert_butyl_pattern = Chem.MolFromSmarts("[#6]-[#6]([#6])([#6])[#6]")
                if mol.HasSubstructMatch(tert_butyl_pattern):
                    tert_butyl_at_depths.add(depth)
                    print(f"Detected tert-butyl group at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if tert-butyl is present at all depths
    all_depths = set(range(max_depth + 1))
    missing_depths = all_depths - tert_butyl_at_depths

    # Return True if tert-butyl is maintained throughout (allowing for some missing depths due to reaction nodes)
    return (
        len(missing_depths) <= max_depth // 2
    )  # Allow for some missing depths due to reaction nodes
