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
    This function detects if the synthetic route involves a late-stage Boc deprotection.
    Late-stage means the deprotection occurs in the final step or penultimate step.
    """
    boc_deprotection_found = False
    depth_of_deprotection = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal boc_deprotection_found, depth_of_deprotection, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check for Boc group in reactants
            reactants_mol = Chem.MolFromSmiles(reactants_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactants_mol and product_mol:
                # Boc group SMARTS pattern
                boc_pattern = Chem.MolFromSmarts("[C;H3][C]([C;H3])([C;H3])[O][C](=[O])[N]")

                # Check if Boc is in reactants but not in product
                if reactants_mol.HasSubstructMatch(
                    boc_pattern
                ) and not product_mol.HasSubstructMatch(boc_pattern):
                    boc_deprotection_found = True
                    depth_of_deprotection = depth
                    print(f"Boc deprotection found at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if deprotection is late-stage (at depth 0 or 1)
    is_late_stage = boc_deprotection_found and depth_of_deprotection <= 1

    if is_late_stage:
        print("Late-stage Boc deprotection strategy detected")

    return is_late_stage
