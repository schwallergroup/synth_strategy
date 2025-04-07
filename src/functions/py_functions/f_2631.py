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
    Detects if the synthesis route involves fragment joining via amide coupling.
    This is a common strategy for connecting molecular fragments.
    """
    amide_coupling_detected = False

    def dfs_traverse(node):
        nonlocal amide_coupling_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for carboxylic acid and amine patterns in reactants and amide in product
            carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=O)[OH]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            amide_pattern = Chem.MolFromSmarts("[C](=O)[NH]")

            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactants_mol and product_mol:
                # Check if carboxylic acid and amine are in reactants and amide is in product
                if (
                    reactants_mol.HasSubstructMatch(carboxylic_acid_pattern)
                    and reactants_mol.HasSubstructMatch(amine_pattern)
                    and product_mol.HasSubstructMatch(amide_pattern)
                ):
                    print("Amide coupling for fragment joining detected")
                    amide_coupling_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_coupling_detected
