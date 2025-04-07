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
    This function detects a linear synthesis with sequential nitrogen acylations
    (amide, carbamate, urea) and nitro reduction.
    """
    # Track key transformations
    acylation_types = []
    nitro_reduction_found = False
    linear_synthesis = True

    def dfs_traverse(node):
        nonlocal acylation_types, nitro_reduction_found, linear_synthesis

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for linear synthesis (only one non-starting material reactant)
            # This is a simplified check - a more robust implementation would track molecule history
            if len(reactants) > 2:
                linear_synthesis = False

            # Create RDKit mol objects
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if all(reactant_mols) and product_mol:
                # Check for nitro reduction
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
                amine_pattern = Chem.MolFromSmarts("[NH2]")

                reactants_with_nitro = any(
                    mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols
                )
                product_has_amine = product_mol.HasSubstructMatch(amine_pattern)

                if (
                    reactants_with_nitro
                    and product_has_amine
                    and not product_mol.HasSubstructMatch(nitro_pattern)
                ):
                    print("Detected nitro reduction to amine")
                    nitro_reduction_found = True

                # Check for acylation products
                amide_pattern = Chem.MolFromSmarts("[#6](=O)[#7]")
                urea_pattern = Chem.MolFromSmarts("[#7][#6](=O)[#7]")
                carbamate_pattern = Chem.MolFromSmarts("[#7][#6](=O)[#8][#6]")

                # Check reactants to see if they already had these patterns
                reactants_with_amide = any(
                    mol.HasSubstructMatch(amide_pattern) for mol in reactant_mols
                )
                reactants_with_urea = any(
                    mol.HasSubstructMatch(urea_pattern) for mol in reactant_mols
                )
                reactants_with_carbamate = any(
                    mol.HasSubstructMatch(carbamate_pattern) for mol in reactant_mols
                )

                # Check for new formations
                if (
                    product_mol.HasSubstructMatch(amide_pattern)
                    and not reactants_with_amide
                ):
                    print("Detected amide formation")
                    acylation_types.append("amide")
                if (
                    product_mol.HasSubstructMatch(urea_pattern)
                    and not reactants_with_urea
                ):
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

    # Check if we have the complete strategy
    unique_acylations = set(acylation_types)
    print(f"Found acylation types: {unique_acylations}")
    print(f"Nitro reduction found: {nitro_reduction_found}")
    print(f"Linear synthesis: {linear_synthesis}")

    return (
        len(unique_acylations) >= 2
        and nitro_reduction_found
        and linear_synthesis
        and "amide" in unique_acylations
    )
