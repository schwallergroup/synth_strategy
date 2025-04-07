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
    This function detects ester hydrolysis to form a carboxylic acid.
    """
    ester_hydrolysis_detected = False

    def dfs_traverse(node):
        nonlocal ester_hydrolysis_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for ester in reactants and carboxylic acid in product
            ester_pattern = Chem.MolFromSmarts("[C;$(C=O)][O][C]")
            acid_pattern = Chem.MolFromSmarts("[C;$(C=O)][OH]")

            reactants_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and acid_pattern:
                has_ester = any(
                    r and r.HasSubstructMatch(ester_pattern)
                    for r in reactants_mols
                    if r
                )
                product_has_acid = product_mol.HasSubstructMatch(acid_pattern)

                if has_ester and product_has_acid:
                    ester_hydrolysis_detected = True
                    print("Ester hydrolysis to carboxylic acid detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return ester_hydrolysis_detected
