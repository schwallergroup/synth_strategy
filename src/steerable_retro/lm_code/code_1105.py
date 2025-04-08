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
    Detects if the synthesis uses amide formation as a key bond-forming step.
    """
    amide_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for carboxylic acid in reactants
            carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

            # Check for amine in reactants
            amine_pattern = Chem.MolFromSmarts("[NX3;H1,H2]")

            # Check for amide in product
            amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")

            has_carboxylic_acid = False
            has_amine = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(carboxylic_acid_pattern):
                    has_carboxylic_acid = True
                if mol and mol.HasSubstructMatch(amine_pattern):
                    has_amine = True

            product_mol = Chem.MolFromSmiles(product_smiles)
            has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

            if has_carboxylic_acid and has_amine and has_amide:
                print(f"Amide formation detected at depth {depth}")
                amide_formation_detected = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return amide_formation_detected
