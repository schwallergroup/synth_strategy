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
    Detects if a nitro group persists through multiple steps of the synthesis.
    """
    nitro_reactions = []

    def dfs_traverse(node, depth=0):
        if node.get("type") == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro group in product
            nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")

            if nitro_pattern:
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                    nitro_reactions.append(depth)
                    print(f"Nitro group found at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if nitro group persists through at least 2 reactions
    persists = len(nitro_reactions) >= 2

    if persists:
        print(f"Nitro group persists through {len(nitro_reactions)} reactions")

    return persists
