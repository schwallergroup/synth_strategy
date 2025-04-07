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
    This function detects if the synthesis involves an amine protection-deprotection sequence,
    specifically looking for Cbz (carbamate) protection.
    """
    protection_depth = None
    deprotection_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal protection_depth, deprotection_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for carbamate protection (amine + Cbz-Cl → N-Cbz)
            product_mol = Chem.MolFromSmiles(product_smiles)
            carbamate_pattern = Chem.MolFromSmarts("[N][C](=[O])[O][C]")

            # Check for amine in reactants
            amine_in_reactants = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[NH2]")
                ):
                    amine_in_reactants = True
                    break

            # Protection: amine → carbamate
            if (
                product_mol
                and product_mol.HasSubstructMatch(carbamate_pattern)
                and amine_in_reactants
            ):
                protection_depth = depth
                print(f"Amine protection detected at depth {depth}")

            # Deprotection: carbamate → amine
            carbamate_in_reactants = False
            for reactant in reactants_smiles:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(carbamate_pattern):
                    carbamate_in_reactants = True
                    break

            product_mol = Chem.MolFromSmiles(product_smiles)
            if (
                product_mol
                and product_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]"))
                and carbamate_in_reactants
            ):
                deprotection_depth = depth
                print(f"Amine deprotection detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if both protection and deprotection were found, and protection happened before deprotection
    has_protection_deprotection = (
        protection_depth is not None
        and deprotection_depth is not None
        and protection_depth > deprotection_depth
    )

    print(f"Amine protection-deprotection sequence: {has_protection_deprotection}")
    print(
        f"Protection depth: {protection_depth}, Deprotection depth: {deprotection_depth}"
    )
    return has_protection_deprotection
