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
    This function detects if the synthesis route involves nitro reduction to amine.
    """
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    amine_pattern = Chem.MolFromSmarts("[NH2]")
    nitro_reduction_detected = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro in reactants and amine in product
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol is not None and any(mol is not None for mol in reactant_mols):
                reactants_have_nitro = any(
                    mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols if mol is not None
                )
                product_has_amine = product_mol.HasSubstructMatch(amine_pattern)

                if reactants_have_nitro and product_has_amine:
                    # This is a simplification - ideally we would check that the amine is at the same position as the nitro
                    nitro_reduction_detected = True
                    print("Nitro reduction detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nitro_reduction_detected
