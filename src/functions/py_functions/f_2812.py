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
    This function detects if the synthesis route contains an alcohol to chloride
    conversion step.
    """
    conversion_found = False

    def dfs_traverse(node):
        nonlocal conversion_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for alcohol to chloride conversion
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol and any(mol for mol in reactant_mols):
                # Check if reactants contain alcohol and product contains chloride
                alcohol_pattern = Chem.MolFromSmarts("[#6][OH]")
                chloride_pattern = Chem.MolFromSmarts("[#6][Cl]")

                has_alcohol = any(
                    mol and mol.HasSubstructMatch(alcohol_pattern)
                    for mol in reactant_mols
                )
                has_chloride = product_mol.HasSubstructMatch(chloride_pattern)

                if has_alcohol and has_chloride:
                    # Additional check to ensure it's a direct conversion
                    # by comparing the molecules more carefully
                    for mol in reactant_mols:
                        if mol and mol.HasSubstructMatch(alcohol_pattern):
                            # This is a simplification - in a real implementation,
                            # you would need more sophisticated matching
                            conversion_found = True
                            print("Alcohol to chloride conversion detected")
                            break

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return conversion_found
