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
    This function detects nitro reduction to amine in the synthesis.
    """
    nitro_reduction_found = False

    def dfs_traverse(node):
        nonlocal nitro_reduction_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            nitro_pattern = Chem.MolFromSmarts("[#6]-[N+](=[O])[O-]")

            # Check for amine group in product
            amine_pattern = Chem.MolFromSmarts("[#6]-[NH2]")

            # Check if reactants have nitro but product has amine instead
            reactants_have_nitro = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(nitro_pattern) for r in reactants if r
            )
            product_has_amine = (
                Chem.MolFromSmiles(product).HasSubstructMatch(amine_pattern) if product else False
            )

            if reactants_have_nitro and product_has_amine:
                nitro_reduction_found = True
                print("Nitro reduction to amine detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitro reduction to amine detected: {nitro_reduction_found}")
    return nitro_reduction_found
