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
    Detects amide formation between a carboxylic acid and an amine.
    """
    found_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_amide_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Patterns for carboxylic acid, amine, and amide
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[OH]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            amide_pattern = Chem.MolFromSmarts("[C](=[O])[NH]")

            # Check for amide formation
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            if (
                any(
                    mol and mol.HasSubstructMatch(carboxylic_acid_pattern)
                    for mol in reactant_mols
                )
                and any(
                    mol and mol.HasSubstructMatch(amine_pattern)
                    for mol in reactant_mols
                )
                and product_mol
                and product_mol.HasSubstructMatch(amide_pattern)
            ):
                found_amide_formation = True
                print("Found amide formation")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Amide formation strategy detected: {found_amide_formation}")
    return found_amide_formation
