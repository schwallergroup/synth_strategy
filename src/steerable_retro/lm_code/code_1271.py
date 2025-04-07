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
    This function detects a synthetic strategy involving nitro group reduction to amine.
    """
    nitro_reduction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_detected

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for nitro reduction pattern
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            reactants_with_nitro = any(r and r.HasSubstructMatch(nitro_pattern) for r in reactants)

            if (
                reactants_with_nitro
                and product
                and product.HasSubstructMatch(amine_pattern)
                and not product.HasSubstructMatch(nitro_pattern)
            ):
                print(f"Detected nitro reduction at depth {depth}")
                nitro_reduction_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Nitro reduction to amine detected: {nitro_reduction_detected}")
    return nitro_reduction_detected
