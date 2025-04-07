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
    Detects a synthesis that includes a nitrile reduction to form an amine.
    """
    has_nitrile_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitrile_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for nitrile reduction
            nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
            amine_pattern = Chem.MolFromSmarts("[C][N;H2]")

            has_nitrile = any(
                r is not None and r.HasSubstructMatch(nitrile_pattern)
                for r in reactants
            )
            has_amine_product = product is not None and product.HasSubstructMatch(
                amine_pattern
            )

            if has_nitrile and has_amine_product:
                has_nitrile_reduction = True
                print("Detected nitrile reduction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_nitrile_reduction
