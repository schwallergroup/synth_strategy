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
    This function detects if the synthetic route contains a Boc protection followed by deprotection.
    """
    boc_protected_intermediates = []
    boc_deprotection_reactions = []

    def dfs_traverse(node):
        nonlocal boc_protected_intermediates, boc_deprotection_reactions

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]
                product_mol = Chem.MolFromSmiles(product)

                # Boc group pattern
                boc_pattern = Chem.MolFromSmarts("[C][C]([C])([C])[O][C](=[O])[N]")

                # Check for Boc deprotection
                has_boc_reactant = any(
                    mol and mol.HasSubstructMatch(boc_pattern) for mol in reactant_mols
                )
                has_boc_product = product_mol and product_mol.HasSubstructMatch(
                    boc_pattern
                )

                if has_boc_reactant and not has_boc_product:
                    boc_deprotection_reactions.append(node)
                    print(f"Found Boc deprotection reaction: {rsmi}")

                # Track Boc-protected intermediates
                if has_boc_product:
                    boc_protected_intermediates.append(product)
                    print(f"Found Boc-protected intermediate: {product}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have both Boc protection and deprotection
    result = (
        len(boc_protected_intermediates) > 0 and len(boc_deprotection_reactions) > 0
    )
    print(f"Boc protection/deprotection sequence detected: {result}")
    return result
