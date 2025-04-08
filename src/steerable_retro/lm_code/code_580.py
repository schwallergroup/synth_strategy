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
    This function detects a synthetic strategy involving late-stage urea formation.
    """
    has_urea_formation = False
    urea_formation_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal has_urea_formation, urea_formation_depth, max_depth

        # Update max depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for urea formation
                urea_pattern = Chem.MolFromSmarts("[#7]-[#6](=[#8])-[#7]")
                product_has_urea = product and product.HasSubstructMatch(urea_pattern)

                if product_has_urea and not any(
                    r and r.HasSubstructMatch(urea_pattern) for r in reactants if r
                ):
                    has_urea_formation = True
                    urea_formation_depth = depth
                    print(f"Found urea formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # For late-stage, the urea formation should be at a low depth (closer to the final product)
    strategy_present = (
        has_urea_formation and urea_formation_depth is not None and urea_formation_depth <= 1
    )

    print(f"Late-stage urea formation strategy detected: {strategy_present}")
    return strategy_present
