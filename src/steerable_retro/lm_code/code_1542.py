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
    This function detects a strategy involving amide formation from carboxylic acid and amine.
    """
    # Initialize tracking variable
    has_amide_formation = False

    def dfs_traverse(node):
        nonlocal has_amide_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation
            carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            amide_pattern = Chem.MolFromSmarts("C(=O)N")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Check if reactants include carboxylic acid and amine
            has_carboxylic_acid = any(
                mol is not None and mol.HasSubstructMatch(carboxylic_acid_pattern)
                for mol in reactant_mols
            )
            has_amine = any(
                mol is not None and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
            )

            # Check if product has amide
            product_has_amide = product_mol is not None and product_mol.HasSubstructMatch(
                amide_pattern
            )

            if has_carboxylic_acid and has_amine and product_has_amide:
                has_amide_formation = True
                print("Detected amide formation reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_amide_formation
