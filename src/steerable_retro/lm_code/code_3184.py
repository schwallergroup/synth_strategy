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
    Detects if the synthesis route involves nitro reduction to amine.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro reduction
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            # Look for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[#6][N+](=O)[O-]")
            has_nitro = any(
                mol is not None and mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols
            )

            # Look for amine in product at same position
            if has_nitro and product_mol is not None:
                amine_pattern = Chem.MolFromSmarts("[#6][NH2]")
                if product_mol.HasSubstructMatch(amine_pattern):
                    has_nitro_reduction = True
                    print("Detected nitro reduction to amine")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return has_nitro_reduction
