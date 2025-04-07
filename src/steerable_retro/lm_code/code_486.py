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
    This function detects if a sulfonamide group is introduced in the second half of the synthesis.
    """
    sulfonamide_formation_depth = None
    max_depth = 0

    # SMARTS patterns
    sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])-[#7]")
    sulfonyl_chloride_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])-[Cl]")

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formation_depth, max_depth

        # Track maximum depth
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            products_part = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
            product = Chem.MolFromSmiles(products_part)

            # Check for sulfonamide formation
            if (
                any(
                    r is not None and r.HasSubstructMatch(sulfonyl_chloride_pattern)
                    for r in reactants
                )
                and product is not None
                and product.HasSubstructMatch(sulfonamide_pattern)
            ):
                sulfonamide_formation_depth = depth
                print(f"Found sulfonamide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if sulfonamide formation occurs in the second half of synthesis
    # (lower depth values correspond to later stages in synthesis)
    if sulfonamide_formation_depth is not None:
        is_late_stage = sulfonamide_formation_depth < (max_depth / 2)
        print(
            f"Sulfonamide formation at depth {sulfonamide_formation_depth}, max depth {max_depth}"
        )
        print(f"Is late stage: {is_late_stage}")
        return is_late_stage

    return False
