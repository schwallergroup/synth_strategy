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
    This function detects if the synthesis includes a reductive amination step
    (aldehyde/ketone + amine → secondary/tertiary amine).
    """
    contains_reductive_amination = False

    def dfs_traverse(node):
        nonlocal contains_reductive_amination

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for aldehyde pattern in reactants
            aldehyde_pattern = Chem.MolFromSmarts("[#6][CH]=O")
            amine_pattern = Chem.MolFromSmarts("[NH2][#6]")

            # Check for secondary amine pattern in product
            sec_amine_pattern = Chem.MolFromSmarts("[#6][CH2][NH][#6]")

            # Check if reactants contain aldehyde and primary amine
            has_aldehyde = False
            has_amine = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(aldehyde_pattern):
                        has_aldehyde = True
                    if mol.HasSubstructMatch(amine_pattern):
                        has_amine = True

            # Check if product contains secondary amine
            product_mol = Chem.MolFromSmiles(product_smiles)
            has_sec_amine = False
            if product_mol and product_mol.HasSubstructMatch(sec_amine_pattern):
                has_sec_amine = True

            # If all conditions are met, this is a reductive amination
            if has_aldehyde and has_amine and has_sec_amine:
                contains_reductive_amination = True
                print("Detected reductive amination strategy")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return contains_reductive_amination
