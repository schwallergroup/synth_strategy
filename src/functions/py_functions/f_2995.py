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
    Detects if the synthesis involves a late-stage nucleophilic aromatic substitution
    on a pyrimidine ring (depth 0 or 1).
    """
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal result

        if node["type"] == "reaction" and depth <= 1:
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if one reactant contains a chloropyrimidine
            pyrimidine_pattern = Chem.MolFromSmarts("[n]1[c][n][c]([Cl,F,Br,I])[n][c]1")
            amine_pattern = Chem.MolFromSmarts("[NH,NH2]")

            has_pyrimidine = False
            has_amine = False

            for reactant in reactants_smiles:
                try:
                    mol = Chem.MolFromSmiles(reactant)
                    if mol and mol.HasSubstructMatch(pyrimidine_pattern):
                        has_pyrimidine = True
                    if mol and mol.HasSubstructMatch(amine_pattern):
                        has_amine = True
                except:
                    continue

            # Check if product has a C-N bond where the chlorine was
            if has_pyrimidine and has_amine:
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(
                        Chem.MolFromSmarts("[n]1[c][n][c]([N])[n][c]1")
                    ):
                        print(
                            "Detected late-stage pyrimidine SNAr reaction at depth",
                            depth,
                        )
                        result = True
                except:
                    pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return result
