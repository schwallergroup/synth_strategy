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
    Detects if the synthesis route involves cyanation via an aryl halide intermediate.
    This is a specific C-C bond formation strategy.
    """
    cyanation_detected = False

    def dfs_traverse(node):
        nonlocal cyanation_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for aryl halide pattern in reactants and aryl nitrile in products
            aryl_iodide_pattern = Chem.MolFromSmarts("c[I]")
            aryl_nitrile_pattern = Chem.MolFromSmarts("cC#N")

            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactants_mol and product_mol:
                # Check if aryl iodide is in reactants and aryl nitrile is in products
                if (
                    reactants_mol.HasSubstructMatch(aryl_iodide_pattern)
                    and product_mol.HasSubstructMatch(aryl_nitrile_pattern)
                    and not reactants_mol.HasSubstructMatch(aryl_nitrile_pattern)
                ):
                    print("Cyanation via aryl halide detected")
                    cyanation_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return cyanation_detected
