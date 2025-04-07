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
    This function detects a synthetic strategy involving late-stage incorporation
    of heterocycles, particularly thiadiazole in the final step.
    """
    depth_of_heterocycle_addition = None
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal depth_of_heterocycle_addition, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for thiadiazole pattern
                thiadiazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#6][#7][#16]1")

                # Check if thiadiazole is in reactants but not in product (addition)
                has_thiadiazole_reactant = any(
                    r.HasSubstructMatch(thiadiazole_pattern) for r in reactants if r
                )
                has_thiadiazole_product = product.HasSubstructMatch(thiadiazole_pattern)

                if has_thiadiazole_reactant and has_thiadiazole_product:
                    if (
                        depth_of_heterocycle_addition is None
                        or depth < depth_of_heterocycle_addition
                    ):
                        depth_of_heterocycle_addition = depth
                        print(f"Detected heterocycle incorporation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it late-stage if it happens in the first half of the synthesis
    # (remember depth 0 is the final product)
    return (
        depth_of_heterocycle_addition is not None
        and depth_of_heterocycle_addition <= max_depth / 2
    )
