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
    This function detects routes with multiple amide coupling reactions.
    """
    amide_coupling_count = 0

    def dfs_traverse(node):
        nonlocal amide_coupling_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for amide coupling: look for amine in reactants and amide in product
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            amine_pattern = Chem.MolFromSmarts("[NH2]")
            carboxylic_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
            acyl_chloride_pattern = Chem.MolFromSmarts("[C](=O)[Cl]")
            amide_pattern = Chem.MolFromSmarts("[C](=O)[NH]")

            # Check if any reactant has amine
            has_amine = any(
                mol.HasSubstructMatch(amine_pattern) for mol in reactants if mol is not None
            )
            # Check if any reactant has carboxylic acid or acyl chloride
            has_carboxylic = any(
                mol.HasSubstructMatch(carboxylic_pattern) for mol in reactants if mol is not None
            )
            has_acyl_chloride = any(
                mol.HasSubstructMatch(acyl_chloride_pattern) for mol in reactants if mol is not None
            )
            # Check if product has amide
            has_amide_product = (
                product.HasSubstructMatch(amide_pattern) if product is not None else False
            )

            # If reactants have amine and (carboxylic acid or acyl chloride) and product has amide, it's likely an amide coupling
            if has_amine and (has_carboxylic or has_acyl_chloride) and has_amide_product:
                amide_coupling_count += 1
                print(f"Detected amide coupling at depth {node.get('depth', 'unknown')}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple amide couplings are detected
    return amide_coupling_count >= 2
