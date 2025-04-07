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
from steerable_retro.utils import check, fuzzy_dict
from steerable_retro.utils.check import Check

fg_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


def main(route):
    """
    Detects a specific sequence of functional group interconversions:
    aldehyde → alcohol → mesylate → nitrile
    """
    # Track transformations we've seen with their connected molecules
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for functional groups in this molecule (for debugging)
            if checker.check_fg("Nitrile", mol_smiles):
                print(f"Found molecule with nitrile at depth {depth}: {mol_smiles}")

            if checker.check_fg("Mesylate", mol_smiles):
                print(f"Found molecule with mesylate at depth {depth}: {mol_smiles}")

            if checker.check_fg("Primary alcohol", mol_smiles) or checker.check_fg(
                "Secondary alcohol", mol_smiles
            ):
                print(f"Found molecule with alcohol at depth {depth}: {mol_smiles}")

            if checker.check_fg("Aldehyde", mol_smiles):
                print(f"Found molecule with aldehyde at depth {depth}: {mol_smiles}")

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                products = rsmi.split(">")[-1].split(".")

                # Check for mesylate to nitrile transformation
                mesylate_reactants = [
                    r for r in reactants if checker.check_fg("Mesylate", r)
                ]
                nitrile_products = [
                    p for p in products if checker.check_fg("Nitrile", p)
                ]

                if mesylate_reactants and nitrile_products:
                    transformations.append(
                        (
                            "mesylate_to_nitrile",
                            depth,
                            mesylate_reactants[0],
                            nitrile_products[0],
                        )
                    )
                    print(
                        f"Found mesylate to nitrile transformation at depth {depth}: {rsmi}"
                    )

                # Check for alcohol to mesylate transformation
                alcohol_reactants = [
                    r
                    for r in reactants
                    if checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                ]
                mesylate_products = [
                    p for p in products if checker.check_fg("Mesylate", p)
                ]

                if alcohol_reactants and mesylate_products:
                    transformations.append(
                        (
                            "alcohol_to_mesylate",
                            depth,
                            alcohol_reactants[0],
                            mesylate_products[0],
                        )
                    )
                    print(
                        f"Found alcohol to mesylate transformation at depth {depth}: {rsmi}"
                    )

                # Check for aldehyde to alcohol transformation
                aldehyde_reactants = [
                    r for r in reactants if checker.check_fg("Aldehyde", r)
                ]
                alcohol_products = [
                    p
                    for p in products
                    if checker.check_fg("Primary alcohol", p)
                    or checker.check_fg("Secondary alcohol", p)
                ]

                if aldehyde_reactants and alcohol_products:
                    transformations.append(
                        (
                            "aldehyde_to_alcohol",
                            depth,
                            aldehyde_reactants[0],
                            alcohol_products[0],
                        )
                    )
                    print(
                        f"Found aldehyde to alcohol transformation at depth {depth}: {rsmi}"
                    )

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth to check sequence
    transformations.sort(key=lambda x: x[1], reverse=True)
    print(f"Sorted transformations: {transformations}")

    # Check if we have all three transformations
    if len(transformations) < 3:
        print("Not all transformations found")
        return False

    # Expected sequence in synthetic order
    expected_sequence = [
        "aldehyde_to_alcohol",
        "alcohol_to_mesylate",
        "mesylate_to_nitrile",
    ]

    # Extract just the transformation types
    transformation_types = [t[0] for t in transformations]
    print(f"Transformation types: {transformation_types}")

    # Check if all expected transformations exist
    for expected in expected_sequence:
        if expected not in transformation_types:
            print(f"Missing transformation: {expected}")
            return False

    # Check if they appear in the correct order
    # In retrosynthetic traversal (higher depth = earlier in synthesis)
    # we expect aldehyde_to_alcohol at highest depth, then alcohol_to_mesylate, then mesylate_to_nitrile

    # Get indices of transformations in our sorted list
    aldehyde_to_alcohol_idx = transformation_types.index("aldehyde_to_alcohol")
    alcohol_to_mesylate_idx = transformation_types.index("alcohol_to_mesylate")
    mesylate_to_nitrile_idx = transformation_types.index("mesylate_to_nitrile")

    # Check correct order by depth
    if not (
        transformations[aldehyde_to_alcohol_idx][1]
        >= transformations[alcohol_to_mesylate_idx][1]
        >= transformations[mesylate_to_nitrile_idx][1]
    ):
        print("Transformations not in correct depth order")
        print(f"Aldehyde→Alcohol depth: {transformations[aldehyde_to_alcohol_idx][1]}")
        print(f"Alcohol→Mesylate depth: {transformations[alcohol_to_mesylate_idx][1]}")
        print(f"Mesylate→Nitrile depth: {transformations[mesylate_to_nitrile_idx][1]}")
        return False

    print("Found complete sequence: aldehyde → alcohol → mesylate → nitrile")
    return True
