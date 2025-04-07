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
    Detects synthesis routes that include nitro group reduction to amine.
    """
    found_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal found_nitro_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro reduction (NO2 → NH2)
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            product_mol = Chem.MolFromSmiles(product_smiles)
            reactant_mols = [
                Chem.MolFromSmiles(r) for r in reactants_smiles if Chem.MolFromSmiles(r)
            ]

            if product_mol and product_mol.HasSubstructMatch(amine_pattern):
                if any(
                    r and r.HasSubstructMatch(nitro_pattern)
                    for r in reactant_mols
                    if r is not None
                ):
                    print("Found nitro reduction")
                    found_nitro_reduction = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_nitro_reduction
