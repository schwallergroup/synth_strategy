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
    This function detects triazole ring formation in the synthetic route.
    """
    triazole_formation_found = False

    def dfs_traverse(node):
        nonlocal triazole_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check if product contains triazole pattern
            product_mol = Chem.MolFromSmiles(product)
            if product_mol and product_mol.HasSubstructMatch(
                Chem.MolFromSmarts("[n]1[n][n][c]2[c]1[n][c][n][c]2")
            ):

                # Check if reactants don't have the complete triazole pattern
                reactants_mol = Chem.MolFromSmiles(reactants)
                if reactants_mol and not reactants_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[n]1[n][n][c]2[c]1[n][c][n][c]2")
                ):
                    print("Triazole formation detected")
                    triazole_formation_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return triazole_formation_found
