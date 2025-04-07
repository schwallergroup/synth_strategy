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
    Detects a sequence of functional group interconversions:
    sulfonic acid → sulfonyl chloride → sulfonamide
    """
    # Track the depths at which each functional group is found
    sulfonic_acid_depth = -1
    sulfonyl_chloride_depth = -1
    sulfonamide_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal sulfonic_acid_depth, sulfonyl_chloride_depth, sulfonamide_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol:
                    # Check for sulfonic acid
                    sulfonic_acid_pattern = Chem.MolFromSmarts(
                        "[c]-[S](=[O])(=[O])-[OH]"
                    )
                    if product_mol.HasSubstructMatch(sulfonic_acid_pattern) and (
                        sulfonic_acid_depth == -1 or depth > sulfonic_acid_depth
                    ):
                        sulfonic_acid_depth = depth
                        print(f"Found sulfonic acid at depth {depth}")

                    # Check for sulfonyl chloride
                    sulfonyl_chloride_pattern = Chem.MolFromSmarts(
                        "[c]-[S](=[O])(=[O])-[Cl]"
                    )
                    if product_mol.HasSubstructMatch(sulfonyl_chloride_pattern) and (
                        sulfonyl_chloride_depth == -1 or depth > sulfonyl_chloride_depth
                    ):
                        sulfonyl_chloride_depth = depth
                        print(f"Found sulfonyl chloride at depth {depth}")

                    # Check for sulfonamide
                    sulfonamide_pattern = Chem.MolFromSmarts("[c]-[S](=[O])(=[O])-[N]")
                    if product_mol.HasSubstructMatch(sulfonamide_pattern) and (
                        sulfonamide_depth == -1 or depth > sulfonamide_depth
                    ):
                        sulfonamide_depth = depth
                        print(f"Found sulfonamide at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if the sequence is found in the correct order (decreasing depth)
    if sulfonic_acid_depth > sulfonyl_chloride_depth > sulfonamide_depth > 0:
        print(
            f"Found sulfonic acid → sulfonyl chloride → sulfonamide sequence: {sulfonic_acid_depth} → {sulfonyl_chloride_depth} → {sulfonamide_depth}"
        )
        return True
    return False
