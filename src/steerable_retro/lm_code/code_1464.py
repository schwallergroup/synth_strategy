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
    This function detects a synthetic strategy where a methylsulfonyl group
    is present throughout the synthesis as a directing/activating group.
    """
    # Initialize tracking variables
    methylsulfonyl_count = 0
    total_reactions = 0

    def dfs_traverse(node):
        nonlocal methylsulfonyl_count, total_reactions

        if node["type"] == "reaction":
            total_reactions += 1

            # Extract product from reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            product = parts[-1]

            # Check if product contains methylsulfonyl group
            product_mol = Chem.MolFromSmiles(product)
            if product_mol:
                methylsulfonyl_pattern = Chem.MolFromSmarts("[CH3][S](=[O])(=[O])[#6]")
                if product_mol.HasSubstructMatch(methylsulfonyl_pattern):
                    methylsulfonyl_count += 1
                    print(f"Detected methylsulfonyl group in reaction product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if methylsulfonyl group is present in all or most reactions
    return (
        total_reactions > 0 and methylsulfonyl_count >= total_reactions * 0.75
    )  # Present in at least 75% of reactions
