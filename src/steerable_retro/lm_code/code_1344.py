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
    Detects if the synthesis route includes an amide formation from carboxylic acid and amine.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for carboxylic acid and amine in reactants
            carboxylic_acid_found = False
            amine_found = False

            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol:
                    carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
                    amine_pattern = Chem.MolFromSmarts("[NH2][c]")

                    if reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                        carboxylic_acid_found = True
                    if reactant_mol.HasSubstructMatch(amine_pattern):
                        amine_found = True

            # Check if product has amide bond
            if carboxylic_acid_found and amine_found:
                product_mol = Chem.MolFromSmiles(product)
                if product_mol:
                    amide_pattern = Chem.MolFromSmarts("[C](=O)[NH][c]")
                    if product_mol.HasSubstructMatch(amide_pattern):
                        amide_formation_found = True
                        print("Amide formation detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_found
