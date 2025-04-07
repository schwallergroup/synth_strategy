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
    Detects late-stage amide coupling (depth 0-1) as the final step in the synthesis.
    """
    has_late_amide_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_amide_coupling

        if node.get("type") == "reaction" and depth <= 1:  # Only check at depths 0-1 (late stage)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carboxylic acid in reactants
                carboxylic_acid_pattern = Chem.MolFromSmarts("O[C](=O)")

                # Check for amine in reactants
                amine_pattern = Chem.MolFromSmarts("[NH]")

                # Check for amide in product
                amide_pattern = Chem.MolFromSmarts("[N]-[C](=O)")

                # Check if reactants contain carboxylic acid and amine, and product contains amide
                has_acid = any(
                    Chem.MolFromSmiles(r)
                    and Chem.MolFromSmiles(r).HasSubstructMatch(carboxylic_acid_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )
                has_amine = any(
                    Chem.MolFromSmiles(r) and Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                    for r in reactants
                    if Chem.MolFromSmiles(r)
                )

                product_mol = Chem.MolFromSmiles(product)
                has_amide = product_mol and product_mol.HasSubstructMatch(amide_pattern)

                if has_acid and has_amine and has_amide:
                    has_late_amide_coupling = True
                    print(f"Found late-stage amide coupling at depth {depth}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_late_amide_coupling
