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
    Detects if the synthesis preserves a benzoxazine core structure throughout.
    """
    benzoxazine_count = 0
    total_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal benzoxazine_count, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for benzoxazine core in product
                product_mol = Chem.MolFromSmiles(product)
                # This is a simplified pattern for benzoxazine-like core
                benzoxazine_pattern = Chem.MolFromSmarts(
                    "[#6]1=[#7][#6]2[#6][#6][#6][#6][#6]2[#6]1"
                )

                if product_mol and product_mol.HasSubstructMatch(benzoxazine_pattern):
                    benzoxazine_count += 1
                    print(f"Benzoxazine core detected in product at depth {depth}")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if benzoxazine core is preserved in most reactions (>70%)
    return total_reactions > 0 and (benzoxazine_count / total_reactions) > 0.7
