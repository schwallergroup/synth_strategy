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
    This function detects if the synthetic route uses a convergent approach
    with 3 or more fragments being combined.
    """
    fragment_count = 0
    visited_nodes = set()

    def count_leaf_nodes(node):
        nonlocal fragment_count, visited_nodes

        # Skip if already visited
        node_id = id(node)
        if node_id in visited_nodes:
            return
        visited_nodes.add(node_id)

        if node["type"] == "mol" and node.get("in_stock", False):
            fragment_count += 1
            print(f"Detected starting material: {node['smiles']}")
            return

        for child in node.get("children", []):
            count_leaf_nodes(child)

    count_leaf_nodes(route)
    is_convergent = fragment_count >= 3
    print(f"Total fragments detected: {fragment_count}")
    return is_convergent
