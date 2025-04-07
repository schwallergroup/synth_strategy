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
    Detects a synthetic strategy involving sequential modifications of a piperazine scaffold.
    """
    # Initialize tracking variables
    piperazine_present = False
    modifications_count = 0

    # SMARTS patterns
    piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")

    def dfs_traverse(node, depth=0):
        nonlocal piperazine_present, modifications_count

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and reactants:
                # Check if piperazine is present
                if product.HasSubstructMatch(piperazine_pattern):
                    piperazine_present = True

                    # Check if this is a modification of the piperazine scaffold
                    reactants_with_piperazine = [
                        r
                        for r in reactants
                        if r and r.HasSubstructMatch(piperazine_pattern)
                    ]

                    if reactants_with_piperazine:
                        # This is a modification of an existing piperazine scaffold
                        modifications_count += 1
                        print(
                            f"Piperazine scaffold modification detected at depth {depth}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if piperazine scaffold is present and modified at least twice
    strategy_present = piperazine_present and modifications_count >= 2
    print(f"Piperazine scaffold modification strategy detected: {strategy_present}")
    return strategy_present
