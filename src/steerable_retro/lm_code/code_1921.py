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
    This function detects if the synthesis includes a nitro group reduction to amine.
    """
    contains_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal contains_nitro_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro pattern in reactants
            nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")

            # Check for amine pattern in product
            amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

            # Check if reactant contains nitro group
            has_nitro = False
            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol and mol.HasSubstructMatch(nitro_pattern):
                    has_nitro = True
                    break

            # Check if product contains amine
            product_mol = Chem.MolFromSmiles(product_smiles)
            has_amine = False
            if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                has_amine = True

            # If both conditions are met, this is a nitro reduction
            if has_nitro and has_amine:
                contains_nitro_reduction = True
                print("Detected nitro to amine reduction strategy")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return contains_nitro_reduction
