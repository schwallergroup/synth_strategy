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
    This function detects if the synthetic route employs multiple alcohol to chloride
    transformations as part of its strategy.
    """
    alcohol_to_chloride_count = 0

    def dfs_traverse(node):
        nonlocal alcohol_to_chloride_count

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to chloride transformation
            alcohol_pattern = Chem.MolFromSmarts("[#6]-[CH2]-[OH]")
            chloride_pattern = Chem.MolFromSmarts("[#6]-[CH2]-[Cl]")

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
            product_mol = Chem.MolFromSmiles(product)

            # Check if any reactant contains alcohol and product contains chloride
            has_alcohol = any(
                mol is not None and mol.HasSubstructMatch(alcohol_pattern) for mol in reactant_mols
            )
            has_chloride = product_mol is not None and product_mol.HasSubstructMatch(
                chloride_pattern
            )

            if has_alcohol and has_chloride:
                print(f"Found alcohol to chloride transformation at reaction with RSMI: {rsmi}")
                alcohol_to_chloride_count += 1

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return alcohol_to_chloride_count >= 2  # At least 2 alcohol to chloride transformations
