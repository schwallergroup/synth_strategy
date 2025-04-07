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
    This function detects a linear synthesis approach with sequential functional group
    transformations (no convergent steps with multiple complex fragments).
    """
    step_count = 0
    is_linear = True
    transformation_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal step_count, is_linear

        if node["type"] == "reaction":
            step_count += 1

            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Count non-trivial reactants (excluding simple reagents)
            complex_reactants = 0
            for r_smiles in reactants_smiles:
                # Skip simple reagents (typically small molecules)
                mol = Chem.MolFromSmiles(r_smiles)
                if mol is not None:
                    atom_count = mol.GetNumAtoms()
                    if atom_count > 5:  # Arbitrary threshold for "complex" molecules
                        complex_reactants += 1

            # If more than 2 complex reactants, it's likely not a linear synthesis
            if complex_reactants > 2:
                is_linear = False
                print(
                    f"Non-linear step detected at depth {depth} with {complex_reactants} complex reactants"
                )

            # Identify transformation type
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product = Chem.MolFromSmiles(product_smiles)

            # Check for specific transformations
            if product is not None:
                # Diaryl ether formation
                diaryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")
                if product.HasSubstructMatch(diaryl_ether_pattern) and not any(
                    r is not None and r.HasSubstructMatch(diaryl_ether_pattern)
                    for r in reactants
                ):
                    transformation_sequence.append(("diaryl_ether_formation", depth))
                    print(f"Detected diaryl ether formation at depth {depth}")

                # Nitrile hydrolysis
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                carboxylic_acid_pattern = Chem.MolFromSmarts("[C](=[O])[O]")
                if any(
                    r is not None and r.HasSubstructMatch(nitrile_pattern)
                    for r in reactants
                ) and product.HasSubstructMatch(carboxylic_acid_pattern):
                    transformation_sequence.append(("nitrile_hydrolysis", depth))
                    print(f"Detected nitrile hydrolysis at depth {depth}")

                # Amide formation
                amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")
                if product.HasSubstructMatch(amide_pattern) and not any(
                    r is not None and r.HasSubstructMatch(amide_pattern)
                    for r in reactants
                ):
                    transformation_sequence.append(("amide_formation", depth))
                    print(f"Detected amide formation at depth {depth}")

                # N-methylation
                n_methyl_pattern = Chem.MolFromSmarts("[N]-[CH3]")
                if product.HasSubstructMatch(n_methyl_pattern) and not any(
                    r is not None and r.HasSubstructMatch(n_methyl_pattern)
                    for r in reactants
                ):
                    transformation_sequence.append(("n_methylation", depth))
                    print(f"Detected N-methylation at depth {depth}")

                # O-methylation
                o_methyl_pattern = Chem.MolFromSmarts("[c]-[O]-[CH3]")
                phenol_pattern = Chem.MolFromSmarts("[c]-[OH]")
                if product.HasSubstructMatch(o_methyl_pattern) and any(
                    r is not None and r.HasSubstructMatch(phenol_pattern)
                    for r in reactants
                ):
                    transformation_sequence.append(("o_methylation", depth))
                    print(f"Detected O-methylation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have a linear synthesis with at least 3 transformations
    has_sufficient_transformations = len(transformation_sequence) >= 3

    # Sort transformations by depth to check sequence
    transformation_sequence.sort(key=lambda x: x[1], reverse=True)

    print(f"Linear synthesis: {is_linear}")
    print(f"Transformation sequence: {transformation_sequence}")

    return is_linear and has_sufficient_transformations
