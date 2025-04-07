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
    This function detects if the synthetic route uses TBDMS protection strategy
    for hydroxyl groups multiple times.
    """
    protection_count = 0

    def dfs_traverse(node):
        nonlocal protection_count

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for TBDMS protection pattern
            # Look for alcohol in reactants
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[#8;H1]")
            # Look for TBDMS ether in product
            tbdms_pattern = Chem.MolFromSmarts(
                "[#8]-[#14]([#6])([#6])[#6]([#6])([#6])[#6]"
            )

            if (
                product
                and any(
                    r for r in reactants if r and r.HasSubstructMatch(alcohol_pattern)
                )
                and product.HasSubstructMatch(tbdms_pattern)
            ):
                protection_count += 1
                print(f"TBDMS protection detected at reaction with SMILES: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple TBDMS protections are found
    result = protection_count >= 2
    print(f"TBDMS protection strategy detected: {result} (count: {protection_count})")
    return result
