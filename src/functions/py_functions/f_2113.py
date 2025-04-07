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
    This function detects if the synthesis route involves hydrolysis of an ester to a carboxylic acid.
    """
    ester_hydrolysis_detected = False

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ester in reactants
            ester_pattern = Chem.MolFromSmarts("C(=O)OC")
            carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")

            # Look for ester in reactants
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(ester_pattern):
                    # Check if the product has a carboxylic acid
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(
                        carboxylic_acid_pattern
                    ):
                        ester_hydrolysis_detected = True
                        print("Ester hydrolysis to carboxylic acid detected")
                        break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return ester_hydrolysis_detected
