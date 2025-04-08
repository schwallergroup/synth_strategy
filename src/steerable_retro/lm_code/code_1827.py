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
    Detects if the synthetic route involves a Boc deprotection in the early stages.
    """
    found_boc_deprotection = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_boc_deprotection, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and reactants:
                # Check for Boc group in reactants but not in product
                boc_pattern = Chem.MolFromSmarts("CC(C)(C)OC(=O)[N]")
                if any(r and r.HasSubstructMatch(boc_pattern) for r in reactants) and (
                    not product.HasSubstructMatch(boc_pattern)
                ):
                    found_boc_deprotection = True
                    print(f"Found Boc deprotection at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Consider it early stage if deprotection occurs in the upper half of the synthesis depth
    early_stage = found_boc_deprotection and (max_depth / 2 < max_depth)
    print(f"Boc deprotection in early stage detected: {early_stage}")
    return early_stage
