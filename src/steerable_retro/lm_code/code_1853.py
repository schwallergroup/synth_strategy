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
    This function detects routes with nitro reduction to amine.
    """
    has_nitro_reduction = False

    def dfs_traverse(node):
        nonlocal has_nitro_reduction

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro reduction: nitro in reactant, amine in product
            for reactant_smiles in reactants_smiles:
                reactant = Chem.MolFromSmiles(reactant_smiles)
                product = Chem.MolFromSmiles(product_smiles)

                if reactant is not None and product is not None:
                    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                    amine_pattern = Chem.MolFromSmarts("[NH2]")

                    if reactant.HasSubstructMatch(nitro_pattern) and product.HasSubstructMatch(
                        amine_pattern
                    ):
                        has_nitro_reduction = True
                        print(f"Detected nitro reduction at depth {node.get('depth', 'unknown')}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_nitro_reduction
