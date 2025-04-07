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
    This function detects a synthetic strategy where a nitro group is preserved
    throughout the synthesis.
    """
    nitro_group_steps = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal nitro_group_steps, total_steps

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            total_steps += 1
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for nitro group in product
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[N+](=[O])[O-]")):
                nitro_group_steps += 1
                print(f"Nitro group detected in step {total_steps}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if nitro group is present in at least half of the steps
    return total_steps > 0 and nitro_group_steps >= total_steps / 2
