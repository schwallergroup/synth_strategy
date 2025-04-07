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
    This function detects a strategy involving nitro reduction to amine,
    followed by nitrogen functionalization ending with sulfonamide formation.
    """
    # Track if we found each transformation
    nitro_reduction_found = False
    amide_formation_found = False
    boc_deprotection_found = False
    sulfonamide_formation_found = False

    # Track the depth of each transformation
    nitro_reduction_depth = -1
    amide_formation_depth = -1
    boc_deprotection_depth = -1
    sulfonamide_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_found, amide_formation_found, boc_deprotection_found, sulfonamide_formation_found
        nonlocal nitro_reduction_depth, amide_formation_depth, boc_deprotection_depth, sulfonamide_formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitro reduction
            nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            reactant_has_nitro = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(nitro_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            product_has_amine = (
                Chem.MolFromSmiles(product).HasSubstructMatch(amine_pattern)
                if Chem.MolFromSmiles(product)
                else False
            )

            if reactant_has_nitro and product_has_amine:
                nitro_reduction_found = True
                nitro_reduction_depth = depth
                print(f"Found nitro reduction at depth {depth}")

            # Check for amide formation
            carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
            amine_pattern = Chem.MolFromSmarts("N")
            amide_pattern = Chem.MolFromSmarts("C(=O)N")

            reactant_has_acid = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(carboxylic_acid_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            reactant_has_amine = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            product_has_amide = (
                Chem.MolFromSmiles(product).HasSubstructMatch(amide_pattern)
                if Chem.MolFromSmiles(product)
                else False
            )

            if reactant_has_acid and reactant_has_amine and product_has_amide:
                amide_formation_found = True
                amide_formation_depth = depth
                print(f"Found amide formation at depth {depth}")

            # Check for Boc deprotection
            boc_pattern = Chem.MolFromSmarts("NC(=O)OC(C)(C)C")
            free_amine_pattern = Chem.MolFromSmarts("N[H]")

            reactant_has_boc = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(boc_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            product_has_free_amine = (
                Chem.MolFromSmiles(product).HasSubstructMatch(free_amine_pattern)
                if Chem.MolFromSmiles(product)
                else False
            )

            if reactant_has_boc and product_has_free_amine:
                boc_deprotection_found = True
                boc_deprotection_depth = depth
                print(f"Found Boc deprotection at depth {depth}")

            # Check for sulfonamide formation
            sulfonyl_chloride_pattern = Chem.MolFromSmarts("S(=O)(=O)Cl")
            amine_pattern = Chem.MolFromSmarts("N")
            sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")

            reactant_has_sulfonyl_chloride = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(sulfonyl_chloride_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            reactant_has_amine = any(
                Chem.MolFromSmiles(r).HasSubstructMatch(amine_pattern)
                for r in reactants
                if Chem.MolFromSmiles(r)
            )
            product_has_sulfonamide = (
                Chem.MolFromSmiles(product).HasSubstructMatch(sulfonamide_pattern)
                if Chem.MolFromSmiles(product)
                else False
            )

            if (
                reactant_has_sulfonyl_chloride
                and reactant_has_amine
                and product_has_sulfonamide
            ):
                sulfonamide_formation_found = True
                sulfonamide_formation_depth = depth
                print(f"Found sulfonamide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # 1. All transformations must be found
    # 2. Sulfonamide formation must be at the lowest depth (latest stage)
    # 3. Nitro reduction must occur before amide formation
    strategy_present = (
        nitro_reduction_found
        and amide_formation_found
        and sulfonamide_formation_found
        and sulfonamide_formation_depth < amide_formation_depth
        and amide_formation_depth < nitro_reduction_depth
    )

    print(f"Strategy present: {strategy_present}")
    return strategy_present
