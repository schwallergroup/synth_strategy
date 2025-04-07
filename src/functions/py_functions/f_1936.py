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
    This function detects if the synthetic route contains heterocycle functionalization,
    particularly focusing on nitrogen-containing heterocycles.
    """
    heterocycle_functionalization_detected = False

    def dfs_traverse(node):
        nonlocal heterocycle_functionalization_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                reactants = Chem.MolFromSmiles(reactants_smiles)
                product = Chem.MolFromSmiles(product_smiles)

                if reactants and product:
                    # Patterns for common nitrogen heterocycles
                    pyrazole_pattern = Chem.MolFromSmarts("[n]1[n][c][c][c]1")
                    triazole_pattern = Chem.MolFromSmarts("[n]1[n][c][n][c]1")

                    # Check for heterocycle functionalization
                    if reactants.HasSubstructMatch(
                        pyrazole_pattern
                    ) or reactants.HasSubstructMatch(triazole_pattern):

                        # Check for nucleophilic substitution on heterocycle
                        if reactants.HasSubstructMatch(
                            Chem.MolFromSmarts("[n][c]([Br,Cl,I,F])")
                        ):
                            if product.HasSubstructMatch(
                                Chem.MolFromSmarts("[n][c]([N])")
                            ):
                                heterocycle_functionalization_detected = True
                                print(f"Heterocycle functionalization detected: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return heterocycle_functionalization_detected
