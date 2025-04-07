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
    Detects a synthetic strategy involving sequential amine protection-deprotection
    with two different protecting groups, followed by sulfonamide formation.
    """
    # Track if we found each key element
    found_trifluoroacetamide_protection = False
    found_trifluoroacetamide_deprotection = False
    found_carbamate_protection = False
    found_carbamate_deprotection = False
    found_sulfonamide_formation = False

    # Track the depth at which each event occurs
    trifluoroacetamide_protection_depth = -1
    trifluoroacetamide_deprotection_depth = -1
    carbamate_protection_depth = -1
    carbamate_deprotection_depth = -1
    sulfonamide_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_trifluoroacetamide_protection, found_trifluoroacetamide_deprotection
        nonlocal found_carbamate_protection, found_carbamate_deprotection
        nonlocal found_sulfonamide_formation
        nonlocal trifluoroacetamide_protection_depth, trifluoroacetamide_deprotection_depth
        nonlocal carbamate_protection_depth, carbamate_deprotection_depth
        nonlocal sulfonamide_formation_depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")

            # Check for trifluoroacetamide protection
            trifluoroacetamide_pattern = Chem.MolFromSmarts("[N]-C(=O)-C([F])([F])([F])")
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol and product_mol.HasSubstructMatch(trifluoroacetamide_pattern):
                # Check if any reactant doesn't have the trifluoroacetamide
                has_reactant_without_trifluoroacetamide = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and not reactant_mol.HasSubstructMatch(
                        trifluoroacetamide_pattern
                    ):
                        has_reactant_without_trifluoroacetamide = True
                        break

                if has_reactant_without_trifluoroacetamide:
                    found_trifluoroacetamide_protection = True
                    trifluoroacetamide_protection_depth = depth
                    print(f"Found trifluoroacetamide protection at depth {depth}")

            # Check for trifluoroacetamide deprotection
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(trifluoroacetamide_pattern):
                    product_mol = Chem.MolFromSmiles(product_part)
                    if product_mol and not product_mol.HasSubstructMatch(
                        trifluoroacetamide_pattern
                    ):
                        found_trifluoroacetamide_deprotection = True
                        trifluoroacetamide_deprotection_depth = depth
                        print(f"Found trifluoroacetamide deprotection at depth {depth}")

            # Check for carbamate protection
            carbamate_pattern = Chem.MolFromSmarts("[N]-C(=O)-O-[C]")
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol and product_mol.HasSubstructMatch(carbamate_pattern):
                # Check if any reactant doesn't have the carbamate
                has_reactant_without_carbamate = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and not reactant_mol.HasSubstructMatch(carbamate_pattern):
                        has_reactant_without_carbamate = True
                        break

                if has_reactant_without_carbamate:
                    found_carbamate_protection = True
                    carbamate_protection_depth = depth
                    print(f"Found carbamate protection at depth {depth}")

            # Check for carbamate deprotection
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(carbamate_pattern):
                    product_mol = Chem.MolFromSmiles(product_part)
                    if product_mol and not product_mol.HasSubstructMatch(carbamate_pattern):
                        found_carbamate_deprotection = True
                        carbamate_deprotection_depth = depth
                        print(f"Found carbamate deprotection at depth {depth}")

            # Check for sulfonamide formation
            sulfonamide_pattern = Chem.MolFromSmarts("[#16](=[#8])(=[#8])-[#7]")
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol and product_mol.HasSubstructMatch(sulfonamide_pattern):
                # Check if any reactant doesn't have the sulfonamide
                has_reactant_without_sulfonamide = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol and not reactant_mol.HasSubstructMatch(sulfonamide_pattern):
                        has_reactant_without_sulfonamide = True
                        break

                if has_reactant_without_sulfonamide:
                    found_sulfonamide_formation = True
                    sulfonamide_formation_depth = depth
                    print(f"Found sulfonamide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    # We need both protection-deprotection cycles and sulfonamide formation
    has_trifluoroacetamide_cycle = (
        found_trifluoroacetamide_protection and found_trifluoroacetamide_deprotection
    )
    has_carbamate_cycle = found_carbamate_protection and found_carbamate_deprotection

    # Check if the sequence is correct (trifluoroacetamide cycle, then carbamate cycle)
    correct_sequence = (
        trifluoroacetamide_protection_depth
        > trifluoroacetamide_deprotection_depth
        > carbamate_protection_depth
        > carbamate_deprotection_depth
    )

    # Check if sulfonamide formation is late-stage (higher depth means earlier in synthesis)
    late_stage_sulfonamide = sulfonamide_formation_depth > trifluoroacetamide_deprotection_depth

    result = (
        has_trifluoroacetamide_cycle
        and has_carbamate_cycle
        and found_sulfonamide_formation
        and correct_sequence
        and late_stage_sulfonamide
    )

    print(f"Sequential amine protection-deprotection with sulfonamide formation: {result}")
    return result
