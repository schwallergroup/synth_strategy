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
    This function detects the reduction of carboxylic acid to alcohol
    in the synthetic route.
    """
    acid_to_alcohol_found = False

    def dfs_traverse(node):
        nonlocal acid_to_alcohol_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid to alcohol transformation
            carboxylic_acid_pattern = Chem.MolFromSmarts("[#6][#6](=[O])[O;H1]")
            alcohol_pattern = Chem.MolFromSmarts("[#6][#6][O;H1]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            has_acid = any(
                mol and mol.HasSubstructMatch(carboxylic_acid_pattern)
                for mol in reactant_mols
                if mol
            )
            has_alcohol = product_mol and product_mol.HasSubstructMatch(alcohol_pattern)

            if has_acid and has_alcohol:
                print("Found carboxylic acid to alcohol transformation")
                acid_to_alcohol_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return acid_to_alcohol_found
