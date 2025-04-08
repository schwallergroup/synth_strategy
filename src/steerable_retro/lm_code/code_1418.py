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
    This function detects a strategy involving sequential sulfonamide formation reactions.
    It looks for two sulfonamide formation reactions (S-N bond formation) in sequence.
    """
    sulfonamide_formation_depths = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a sulfonamide formation reaction
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Look for sulfonyl chloride pattern in reactants
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("S(=O)(=O)Cl")
            has_sulfonyl_chloride = any(
                mol.HasSubstructMatch(sulfonyl_chloride_pattern)
                for mol in reactant_mols
                if mol is not None
            )

            # Look for sulfonamide pattern in product
            sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)[NH]")
            has_sulfonamide_product = product_mol is not None and product_mol.HasSubstructMatch(
                sulfonamide_pattern
            )

            # If both patterns are found, this is likely a sulfonamide formation
            if has_sulfonyl_chloride and has_sulfonamide_product:
                depth = node.get("depth", 0)
                sulfonamide_formation_depths.append(depth)
                print(f"Found sulfonamide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have sequential sulfonamide formations
    # Sort depths to find consecutive reactions
    sulfonamide_formation_depths.sort()

    # Check if we have at least 2 sulfonamide formations
    if len(sulfonamide_formation_depths) >= 2:
        # Check if any two formations are at consecutive or nearby depths
        for i in range(len(sulfonamide_formation_depths) - 1):
            if abs(sulfonamide_formation_depths[i] - sulfonamide_formation_depths[i + 1]) <= 2:
                print("Found sequential sulfonamide formations")
                return True

    return False
