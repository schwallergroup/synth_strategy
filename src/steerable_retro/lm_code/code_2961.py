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
    This function detects the retention of a bromoaryl group throughout the synthesis
    (potential handle for further functionalization).
    """
    # Initialize tracking variables
    bromoaryl_count = 0
    total_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal bromoaryl_count, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1
            # Extract product
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]
            product = Chem.MolFromSmiles(product_smiles)

            # Bromoaryl pattern
            bromoaryl_pattern = Chem.MolFromSmarts("c[Br]")

            # Check if product has bromoaryl group
            if product is not None and product.HasSubstructMatch(bromoaryl_pattern):
                bromoaryl_count += 1
                print(f"Detected bromoaryl group at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if bromoaryl group is retained throughout most/all of the synthesis
    return bromoaryl_count > 0 and bromoaryl_count == total_reactions
