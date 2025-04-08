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
    Detects synthesis routes with multiple esterification reactions
    (COOH → COOMe transformations).
    """
    esterification_count = 0

    def dfs_traverse(node):
        nonlocal esterification_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for esterification (COOH → COOMe)
            carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
            methyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")

            product_mol = Chem.MolFromSmiles(product_smiles)
            reactant_mols = [
                Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
            ]

            # Check if any reactant has carboxylic acid and product has methyl ester
            if product_mol and product_mol.HasSubstructMatch(methyl_ester_pattern):
                if any(
                    r and r.HasSubstructMatch(carboxylic_acid_pattern)
                    for r in reactant_mols
                    if r is not None
                ):
                    print(f"Found esterification at depth: {node.get('depth', 'unknown')}")
                    esterification_count += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return esterification_count >= 2
