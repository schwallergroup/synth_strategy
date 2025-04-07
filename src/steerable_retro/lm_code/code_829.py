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
    This function detects a synthetic strategy involving conversion of aryl hydroxyl
    to aryl chloride.
    """
    has_hydroxyl_to_chloride = False

    def dfs_traverse(node, depth=0):
        nonlocal has_hydroxyl_to_chloride

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aryl hydroxyl in reactants and aryl chloride in product
            aryl_hydroxyl_pattern = Chem.MolFromSmarts("[c][OH]")
            aryl_chloride_pattern = Chem.MolFromSmarts("[c][Cl]")

            has_aryl_hydroxyl = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(aryl_hydroxyl_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            product_mol = Chem.MolFromSmiles(product)
            has_aryl_chloride = product_mol and product_mol.HasSubstructMatch(aryl_chloride_pattern)

            if has_aryl_hydroxyl and has_aryl_chloride:
                print(f"Detected hydroxyl to chloride conversion at depth {depth}")
                has_hydroxyl_to_chloride = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_hydroxyl_to_chloride
