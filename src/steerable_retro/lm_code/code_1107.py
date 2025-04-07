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
    Detects if the synthesis uses Suzuki coupling for aryl-aryl bond formation.
    """
    suzuki_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal suzuki_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for boronic acid in reactants
            boronic_acid_pattern = Chem.MolFromSmarts("[c][B]([OH])[OH]")

            # Check for aryl halide in reactants
            aryl_halide_pattern = Chem.MolFromSmarts("[c][Br,I,Cl]")

            # Check for biaryl in product
            biaryl_pattern = Chem.MolFromSmarts("[c]-[c]")

            has_boronic_acid = False
            has_aryl_halide = False

            for reactant in reactants_smiles:
                mol = Chem.MolFromSmiles(reactant)
                if mol:
                    if mol.HasSubstructMatch(boronic_acid_pattern):
                        has_boronic_acid = True
                    if mol.HasSubstructMatch(aryl_halide_pattern):
                        has_aryl_halide = True

            product_mol = Chem.MolFromSmiles(product_smiles)
            has_biaryl = product_mol and product_mol.HasSubstructMatch(biaryl_pattern)

            if has_boronic_acid and has_aryl_halide and has_biaryl:
                print(f"Suzuki coupling detected at depth {depth}")
                suzuki_detected = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return suzuki_detected
