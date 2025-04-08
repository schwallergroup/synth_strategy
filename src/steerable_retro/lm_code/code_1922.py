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
    This function detects if the synthesis includes a phenol alkylation step
    (phenol + alkyl halide → aryl ether).
    """
    contains_phenol_alkylation = False

    def dfs_traverse(node):
        nonlocal contains_phenol_alkylation

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for phenol pattern in reactants
            phenol_pattern = Chem.MolFromSmarts("c[OH]")
            alkyl_halide_pattern = Chem.MolFromSmarts("[#6][Cl,Br,I]")

            # Check for aryl ether pattern in product
            aryl_ether_pattern = Chem.MolFromSmarts("c[O][#6]")

            # Check if reactants contain phenol and alkyl halide
            has_phenol = False
            has_alkyl_halide = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(phenol_pattern):
                        has_phenol = True
                    if mol.HasSubstructMatch(alkyl_halide_pattern):
                        has_alkyl_halide = True

            # Check if product contains aryl ether
            product_mol = Chem.MolFromSmiles(product_smiles)
            has_aryl_ether = False
            if product_mol and product_mol.HasSubstructMatch(aryl_ether_pattern):
                has_aryl_ether = True

            # If all conditions are met, this is a phenol alkylation
            if has_phenol and has_alkyl_halide and has_aryl_ether:
                contains_phenol_alkylation = True
                print("Detected phenol alkylation strategy")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return contains_phenol_alkylation
