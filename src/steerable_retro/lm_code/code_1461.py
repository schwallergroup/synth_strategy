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
    Detects a synthetic strategy with early biaryl formation (Suzuki coupling)
    and late-stage amide formation.
    """
    # Initialize tracking variables
    has_biaryl_formation = False
    has_amide_formation = False
    biaryl_depth = -1
    amide_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal has_biaryl_formation, has_amide_formation
        nonlocal biaryl_depth, amide_formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product_mol and len(reactant_mols) >= 1:
                # Check for Suzuki coupling (biaryl formation)
                if len(reactant_mols) >= 2:
                    # Look for boronic acid/ester and aryl halide patterns
                    boronic_acid_pattern = Chem.MolFromSmarts("[c]-[B]")
                    aryl_halide_pattern = Chem.MolFromSmarts("[c]-[Cl,Br,I,F]")

                    has_boronic_acid = any(
                        mol.HasSubstructMatch(boronic_acid_pattern) for mol in reactant_mols if mol
                    )
                    has_aryl_halide = any(
                        mol.HasSubstructMatch(aryl_halide_pattern) for mol in reactant_mols if mol
                    )

                    # Check if product has new biaryl bond
                    biaryl_pattern = Chem.MolFromSmarts("c-c")

                    if (
                        has_boronic_acid
                        and has_aryl_halide
                        and product_mol.HasSubstructMatch(biaryl_pattern)
                    ):
                        has_biaryl_formation = True
                        biaryl_depth = depth
                        print(f"Detected biaryl formation at depth {depth}")

                # Check for amide formation
                carboxylic_acid_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[#8]")
                amine_pattern = Chem.MolFromSmarts("[#6]-[#7](-[#6])-[#1]")
                amide_pattern = Chem.MolFromSmarts("[#6]-[#6](=[O])-[#7](-[#6])-[#6]")

                has_acid_reactant = any(
                    mol.HasSubstructMatch(carboxylic_acid_pattern) for mol in reactant_mols if mol
                )
                has_amine_reactant = any(
                    mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols if mol
                )
                has_amide_product = product_mol and product_mol.HasSubstructMatch(amide_pattern)

                if (has_acid_reactant or has_amine_reactant) and has_amide_product:
                    has_amide_formation = True
                    amide_formation_depth = depth
                    print(f"Detected amide formation at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if the strategy is present - early biaryl (higher depth) and late amide (lower depth)
    strategy_present = (
        has_biaryl_formation
        and has_amide_formation
        and biaryl_depth > amide_formation_depth
        and amide_formation_depth <= 1  # Ensure amide formation is late-stage
    )

    if strategy_present:
        print("Detected strategy: Early biaryl formation with late-stage amide coupling")

    return strategy_present
