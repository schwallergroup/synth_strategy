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
    This function detects routes with acyl chloride formation.
    """
    has_acyl_chloride_formation = False

    def dfs_traverse(node):
        nonlocal has_acyl_chloride_formation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for acyl chloride formation: carboxylic acid in reactants, acyl chloride in product
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product is not None:
                carboxylic_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                acyl_chloride_pattern = Chem.MolFromSmarts("[C](=O)[Cl]")

                # Check if any reactant has carboxylic acid
                has_carboxylic = any(
                    mol.HasSubstructMatch(carboxylic_pattern)
                    for mol in reactants
                    if mol is not None
                )
                # Check if product has acyl chloride
                has_acyl_chloride = product.HasSubstructMatch(acyl_chloride_pattern)

                if has_carboxylic and has_acyl_chloride:
                    has_acyl_chloride_formation = True
                    print(
                        f"Detected acyl chloride formation at depth {node.get('depth', 'unknown')}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_acyl_chloride_formation
