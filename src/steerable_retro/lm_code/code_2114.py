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
    This function detects if the synthesis route involves formation of an amide from a carboxylic acid.
    """
    amide_formation_detected = False

    def dfs_traverse(node):
        nonlocal amide_formation_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for carboxylic acid in reactants
            carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
            amide_pattern = Chem.MolFromSmarts("C(=O)N")

            # Look for carboxylic acid in reactants
            acid_in_reactants = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(carboxylic_acid_pattern):
                    acid_in_reactants = True
                    break

            # Check if product has an amide
            product_mol = Chem.MolFromSmiles(product_smiles)
            if acid_in_reactants and product_mol and product_mol.HasSubstructMatch(amide_pattern):
                amide_formation_detected = True
                print("Amide formation from carboxylic acid detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return amide_formation_detected
