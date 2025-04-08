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
    This function detects if the synthetic route includes aldehyde reduction
    to alcohol as a key transformation.
    """
    reduction_found = False

    def dfs_traverse(node):
        nonlocal reduction_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for aldehyde reduction pattern
            aldehyde_pattern = Chem.MolFromSmarts("[#6]=O")
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")

            if (
                product
                and any(r for r in reactants if r and r.HasSubstructMatch(aldehyde_pattern))
                and product.HasSubstructMatch(alcohol_pattern)
            ):
                reduction_found = True
                print(f"Aldehyde reduction detected in reaction: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Aldehyde reduction strategy detected: {reduction_found}")
    return reduction_found
