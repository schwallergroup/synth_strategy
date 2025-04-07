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
    This function detects a strategy involving azide reduction to primary amine.
    """
    has_azide_reduction = False

    def dfs_traverse(node):
        nonlocal has_azide_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an azide reduction reaction
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Look for azide pattern in reactants
            azide_pattern = Chem.MolFromSmarts("[N-]=[N+]=N")
            has_azide = any(
                mol is not None and mol.HasSubstructMatch(azide_pattern)
                for mol in reactant_mols
            )

            # Look for primary amine pattern in product
            amine_pattern = Chem.MolFromSmarts("[NH2]C")
            has_amine_product = (
                product_mol is not None and product_mol.HasSubstructMatch(amine_pattern)
            )

            # If azide is in reactants and amine is in product, this is likely an azide reduction
            if has_azide and has_amine_product:
                has_azide_reduction = True
                print("Found azide reduction to amine")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_azide_reduction
