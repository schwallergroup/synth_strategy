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
    This function detects a sequence of different nitrogen acylation reactions
    (amide, carbamate, urea) in the synthesis route.
    """
    acylation_types = []

    def dfs_traverse(node):
        nonlocal acylation_types

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Create RDKit mol objects
            product_mol = Chem.MolFromSmiles(product) if product else None

            if product_mol:
                # Check for newly formed acylation products
                amide_pattern = Chem.MolFromSmarts("[#6](=O)[#7]")
                urea_pattern = Chem.MolFromSmarts("[#7][#6](=O)[#7]")
                carbamate_pattern = Chem.MolFromSmarts("[#7][#6](=O)[#8][#6]")

                # Check reactants to see if they already had these patterns
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                reactants_with_amide = any(
                    mol and mol.HasSubstructMatch(amide_pattern) for mol in reactant_mols
                )
                reactants_with_urea = any(
                    mol and mol.HasSubstructMatch(urea_pattern) for mol in reactant_mols
                )
                reactants_with_carbamate = any(
                    mol and mol.HasSubstructMatch(carbamate_pattern) for mol in reactant_mols
                )

                # Check for new formations
                if product_mol.HasSubstructMatch(amide_pattern) and not reactants_with_amide:
                    print("Detected amide formation")
                    acylation_types.append("amide")
                if product_mol.HasSubstructMatch(urea_pattern) and not reactants_with_urea:
                    print("Detected urea formation")
                    acylation_types.append("urea")
                if (
                    product_mol.HasSubstructMatch(carbamate_pattern)
                    and not reactants_with_carbamate
                ):
                    print("Detected carbamate formation")
                    acylation_types.append("carbamate")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have at least 2 different types of acylations
    unique_acylations = set(acylation_types)
    print(f"Found acylation types: {unique_acylations}")
    return len(unique_acylations) >= 2
