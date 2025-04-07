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
    This function detects sulfonamide formation from sulfonyl chloride and amine.
    """
    sulfonamide_formation_found = False

    def dfs_traverse(node):
        nonlocal sulfonamide_formation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for sulfonyl chloride in reactants
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("S(=O)(=O)Cl")

            # Check for sulfonamide in product
            sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")

            # Check if reactants have sulfonyl chloride and product has sulfonamide
            reactants_have_sulfonyl_chloride = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(sulfonyl_chloride_pattern)
                for r in reactants
                if r
            )
            product_has_sulfonamide = (
                Chem.MolFromSmiles(product).HasSubstructMatch(sulfonamide_pattern)
                if product
                else False
            )

            if reactants_have_sulfonyl_chloride and product_has_sulfonamide:
                sulfonamide_formation_found = True
                print("Sulfonamide formation detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Sulfonamide formation detected: {sulfonamide_formation_found}")
    return sulfonamide_formation_found
