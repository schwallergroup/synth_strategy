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
    This function detects if the synthesis includes an acetal protection step
    (aldehyde + alcohol → acetal).
    """
    contains_acetal_protection = False

    def dfs_traverse(node):
        nonlocal contains_acetal_protection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aldehyde pattern in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[#6][CH]=O")
            cyclic_ether_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#8]1")

            # Check for acetal pattern in product
            acetal_pattern = Chem.MolFromSmarts("[#6][#8][#6]([#8][#6])[#6]")

            # Check if reactants contain aldehyde and cyclic ether
            has_aldehyde = False
            has_cyclic_ether = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(aldehyde_pattern):
                        has_aldehyde = True
                    if mol.HasSubstructMatch(cyclic_ether_pattern):
                        has_cyclic_ether = True

            # Check if product contains acetal
            product_mol = Chem.MolFromSmiles(product_smiles)
            has_acetal = False
            if product_mol and product_mol.HasSubstructMatch(acetal_pattern):
                has_acetal = True

            # If all conditions are met, this is an acetal protection
            if has_aldehyde and has_cyclic_ether and has_acetal:
                contains_acetal_protection = True
                print("Detected acetal protection strategy")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return contains_acetal_protection
