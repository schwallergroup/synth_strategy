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
    Detects a synthetic strategy involving late-stage nucleophilic substitution
    of a benzyl halide with a cyclic amine, following earlier nitro reduction and amide formation.
    """
    # Track if we found each key transformation
    found_nitro_reduction = False
    found_amide_formation = False
    found_benzyl_halide_substitution = False

    # Track if the benzyl halide substitution is in the late stage (low depth)
    benzyl_halide_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal found_nitro_reduction, found_amide_formation, found_benzyl_halide_substitution, benzyl_halide_depth

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if not product_mol or not all(reactant_mols):
                print("Warning: Could not parse some molecules in reaction")
                return

            # Check for nitro reduction
            nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            if any(
                mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols
            ) and product_mol.HasSubstructMatch(amine_pattern):
                found_nitro_reduction = True
                print(f"Found nitro reduction at depth {depth}")

            # Check for amide formation
            amine_reactant = any(
                mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
            )
            amide_pattern = Chem.MolFromSmarts("[NH]C(=O)")

            if amine_reactant and product_mol.HasSubstructMatch(amide_pattern):
                found_amide_formation = True
                print(f"Found amide formation at depth {depth}")

            # Check for benzyl halide substitution with cyclic amine
            benzyl_halide_pattern = Chem.MolFromSmarts("[c][CH2][Cl,Br,I]")
            cyclic_amine_pattern = Chem.MolFromSmarts(
                "[N]1[C][C][C][C][C]1"
            )  # Piperidine pattern
            benzyl_amine_pattern = Chem.MolFromSmarts("[c][CH2][N]")

            benzyl_halide_present = any(
                mol.HasSubstructMatch(benzyl_halide_pattern) for mol in reactant_mols
            )
            cyclic_amine_present = any(
                mol.HasSubstructMatch(cyclic_amine_pattern) for mol in reactant_mols
            )
            benzyl_amine_in_product = product_mol.HasSubstructMatch(
                benzyl_amine_pattern
            )

            if (
                benzyl_halide_present
                and cyclic_amine_present
                and benzyl_amine_in_product
            ):
                found_benzyl_halide_substitution = True
                benzyl_halide_depth = depth
                print(f"Found benzyl halide substitution at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found all required transformations and if benzyl halide substitution is late-stage
    strategy_present = (
        found_nitro_reduction
        and found_amide_formation
        and found_benzyl_halide_substitution
        and benzyl_halide_depth is not None
        and benzyl_halide_depth <= 1
    )  # Low depth means late-stage

    print(f"Strategy detection result: {strategy_present}")
    return strategy_present
