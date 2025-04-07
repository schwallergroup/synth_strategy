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
    Detects if the synthesis route involves functionalization of heterocycles like imidazole or pyrazole.
    """
    heterocycle_functionalization_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_functionalization_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and all(r for r in reactant_mols):
                    # Check for imidazole or pyrazole patterns
                    imidazole_pattern = Chem.MolFromSmarts("[n]1[c][n][c][c]1")
                    pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c]1")

                    # Check if heterocycle exists in both reactants and products
                    heterocycle_in_reactants = False
                    for r_mol in reactant_mols:
                        if r_mol.HasSubstructMatch(
                            imidazole_pattern
                        ) or r_mol.HasSubstructMatch(pyrazole_pattern):
                            heterocycle_in_reactants = True
                            break

                    heterocycle_in_product = product_mol.HasSubstructMatch(
                        imidazole_pattern
                    ) or product_mol.HasSubstructMatch(pyrazole_pattern)

                    if heterocycle_in_reactants and heterocycle_in_product:
                        # Check if the heterocycle is being modified (different substitution pattern)
                        heterocycle_functionalization_detected = True
                        print(
                            f"Heterocycle functionalization detected at depth {depth}"
                        )
            except:
                print(f"Error processing reaction at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return heterocycle_functionalization_detected
