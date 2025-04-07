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
    Detects a synthetic strategy involving the formation of a sulfonamide linkage
    between two aromatic fragments.
    """
    sulfonamide_formation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formation_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if any reactant has a sulfonyl chloride group
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("S(=O)(=O)Cl")
            has_sulfonyl_chloride = False

            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(sulfonyl_chloride_pattern):
                        has_sulfonyl_chloride = True
                        break
                except:
                    continue

            # Check if product has a sulfonamide group
            sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")
            has_sulfonamide_product = False

            try:
                mol = Chem.MolFromSmiles(product_smiles)
                if mol and mol.HasSubstructMatch(sulfonamide_pattern):
                    has_sulfonamide_product = True
            except:
                pass

            # If both conditions are met, we have a sulfonamide formation
            if has_sulfonyl_chloride and has_sulfonamide_product:
                sulfonamide_formation_found = True
                print(
                    "Found sulfonamide formation: sulfonyl chloride reacted with amine"
                )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return sulfonamide_formation_found
