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
    This function detects ester reduction to alcohol in the synthesis.
    """
    ester_reduction_found = False

    def dfs_traverse(node):
        nonlocal ester_reduction_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ester reduction
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Ester pattern
            ester_pattern = Chem.MolFromSmarts("[#6]-C(=O)O[#6]")

            # Alcohol pattern
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[OH]")

            # Check if ester is in reactants and alcohol is in product
            if product_mol and product_mol.HasSubstructMatch(alcohol_pattern):
                reactants_have_ester = any(
                    mol and mol.HasSubstructMatch(ester_pattern)
                    for mol in reactant_mols
                    if mol
                )
                if reactants_have_ester:
                    # Make sure the alcohol is new (not present in reactants)
                    reactants_have_alcohol = any(
                        mol and mol.HasSubstructMatch(alcohol_pattern)
                        for mol in reactant_mols
                        if mol
                    )
                    if not reactants_have_alcohol:
                        print("Ester reduction to alcohol detected")
                        ester_reduction_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return ester_reduction_found
