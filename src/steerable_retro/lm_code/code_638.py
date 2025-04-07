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
    This function detects a synthetic strategy involving protection/deprotection sequences,
    specifically focusing on Boc protection of amines.
    """
    has_protection = False

    def dfs_traverse(node):
        nonlocal has_protection

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            if product and all(r for r in reactants):
                # Check for Boc protection
                boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
                boc_anhydride_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)OC(=O)O")

                has_boc_product = product.HasSubstructMatch(boc_pattern)
                has_boc_reagent = any(
                    r.HasSubstructMatch(boc_anhydride_pattern) for r in reactants if r
                )

                if has_boc_product and has_boc_reagent:
                    has_protection = True
                    print("Detected Boc protection")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_protection
