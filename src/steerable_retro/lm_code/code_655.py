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
    This function detects nitro to amine reduction as part of the synthetic strategy.
    """
    reduction_found = False

    def dfs_traverse(node):
        nonlocal reduction_found

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for nitro to amine reduction
            nitro_pattern = Chem.MolFromSmarts("[c]-[N+](=[O])-[O-]")
            amine_pattern = Chem.MolFromSmarts("[c]-[NH2]")

            reactants_mol = Chem.MolFromSmiles(reactants)
            product_mol = Chem.MolFromSmiles(product)

            if reactants_mol and product_mol:
                if reactants_mol.HasSubstructMatch(nitro_pattern) and product_mol.HasSubstructMatch(
                    amine_pattern
                ):
                    reduction_found = True
                    print("Nitro to amine reduction found")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return reduction_found
