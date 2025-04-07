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
    This function detects the use of sulfonamide as a protection strategy for amines.
    """
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            product_mol = Chem.MolFromSmiles(product_smiles)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]

            # Check for sulfonamide formation (protection)
            sulfonamide_pattern = Chem.MolFromSmarts("[NH]-[S](=O)(=O)[#6]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")
            sulfonyl_pattern = Chem.MolFromSmarts("[#6]-[S](=O)(=O)[Cl]")

            if (
                product_mol is not None
                and product_mol.HasSubstructMatch(sulfonamide_pattern)
                and any(
                    r is not None and r.HasSubstructMatch(amine_pattern)
                    for r in reactant_mols
                )
                and any(
                    r is not None and r.HasSubstructMatch(sulfonyl_pattern)
                    for r in reactant_mols
                )
            ):
                protection_found = True
                print("Sulfonamide protection detected")

            # Check for sulfonamide cleavage (deprotection)
            if (
                any(
                    r is not None and r.HasSubstructMatch(sulfonamide_pattern)
                    for r in reactant_mols
                )
                and product_mol is not None
                and product_mol.HasSubstructMatch(amine_pattern)
            ):
                deprotection_found = True
                print("Sulfonamide deprotection detected")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Protection strategy involves both protection and deprotection
    strategy_detected = protection_found or deprotection_found

    if strategy_detected:
        print("Protection/deprotection strategy with sulfonamide detected")

    return strategy_detected
