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
    This function detects routes with esterification (carboxylic acid to ester).
    """
    has_esterification = False

    def dfs_traverse(node):
        nonlocal has_esterification

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for esterification: carboxylic acid in reactants, ester in product
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product is not None:
                carboxylic_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                alcohol_pattern = Chem.MolFromSmarts("[C][OH]")
                ester_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")

                # Check if any reactant has carboxylic acid
                has_carboxylic = any(
                    mol.HasSubstructMatch(carboxylic_pattern)
                    for mol in reactants
                    if mol is not None
                )
                # Check if any reactant has alcohol
                has_alcohol = any(
                    mol.HasSubstructMatch(alcohol_pattern) for mol in reactants if mol is not None
                )
                # Check if product has ester
                has_ester_product = product.HasSubstructMatch(ester_pattern)

                if has_carboxylic and has_alcohol and has_ester_product:
                    has_esterification = True
                    print(f"Detected esterification at depth {node.get('depth', 'unknown')}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_esterification
