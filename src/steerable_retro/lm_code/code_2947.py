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
    This function detects if the synthetic route employs carboxylic acid reduction
    as part of its strategy.
    """
    carboxylic_acid_reduction_found = False

    def dfs_traverse(node):
        nonlocal carboxylic_acid_reduction_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid to alcohol transformation
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]-C(=O)[OH]")
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[CH2]-[OH]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Check if any reactant contains carboxylic acid and product contains alcohol
            has_carboxylic_acid = any(
                mol is not None and mol.HasSubstructMatch(carboxylic_acid_pattern)
                for mol in reactant_mols
            )
            has_alcohol = product_mol is not None and product_mol.HasSubstructMatch(alcohol_pattern)

            if has_carboxylic_acid and has_alcohol:
                print(f"Found carboxylic acid reduction at reaction with RSMI: {rsmi}")
                carboxylic_acid_reduction_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return carboxylic_acid_reduction_found
