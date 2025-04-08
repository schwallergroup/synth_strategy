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

root_data = "/home/andres/Documents/steerable_retro/data"

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
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
    This function detects if the synthesis involves a sequence of
    vinyl → aldehyde → heterocycle → amine transformations.
    """
    # Track the sequence of functional groups with their depth
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    print(f"No reaction SMILES found at depth {depth}")
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")
                print(f"Product has vinyl: {checker.check_fg('Vinyl', product_smiles)}")
                print(f"Product has aldehyde: {checker.check_fg('Aldehyde', product_smiles)}")

                # Check for heterocycles
                has_isoxazole = checker.check_ring("isoxazole", product_smiles)
                nitrogen_heterocycles = [
                    "oxazole",
                    "thiazole",
                    "pyrazole",
                    "imidazole",
                    "triazole",
                    "tetrazole",
                ]
                has_other_heterocycle = any(
                    checker.check_ring(ring, product_smiles) for ring in nitrogen_heterocycles
                )

                print(f"Product has isoxazole: {has_isoxazole}")
                print(f"Product has other heterocycle: {has_other_heterocycle}")
                print(
                    f"Product has primary amine: {checker.check_fg('Primary amine', product_smiles)}"
                )

                # 1. Vinyl introduction
                if checker.check_fg("Vinyl", product_smiles):
                    # Check for common vinyl-forming reactions
                    if (
                        checker.check_reaction("Wittig", rsmi)
                        or checker.check_reaction("Wittig reaction with triphenylphosphorane", rsmi)
                        or checker.check_reaction("Wittig with Phosphonium", rsmi)
                        or checker.check_reaction("Julia Olefination", rsmi)
                        or checker.check_reaction("Heck terminal vinyl", rsmi)
                        or checker.check_reaction("Stille reaction_vinyl", rsmi)
                    ):
                        print(
                            f"Vinyl introduction detected with specific reaction at depth {depth}"
                        )
                    else:
                        print(f"Vinyl detected at depth {depth}")

                    # Add to transformations regardless of specific reaction
                    transformations.append(("vinyl", depth))

                # 2. Aldehyde formation
                if checker.check_fg("Aldehyde", product_smiles):
                    if any(checker.check_fg("Vinyl", r) for r in reactants_smiles):
                        if checker.check_reaction("Alkene oxidation to aldehyde", rsmi):
                            print(
                                f"Aldehyde formation from vinyl detected with specific reaction at depth {depth}"
                            )
                        else:
                            print(f"Aldehyde formation detected at depth {depth}")
                    else:
                        print(f"Aldehyde detected at depth {depth}")

                    # Add to transformations regardless of source
                    transformations.append(("aldehyde", depth))

                # 3. Heterocycle formation
                if has_isoxazole:
                    print(f"Isoxazole detected at depth {depth}")
                    transformations.append(("heterocycle", depth))
                elif has_other_heterocycle:
                    for ring in nitrogen_heterocycles:
                        if checker.check_ring(ring, product_smiles):
                            print(f"{ring.capitalize()} detected at depth {depth}")
                            transformations.append(("heterocycle", depth))
                            break

                # 4. Amine formation
                if checker.check_fg("Primary amine", product_smiles):
                    # Check for nitro reduction
                    if any(checker.check_fg("Nitro group", r) for r in reactants_smiles):
                        if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                            print(f"Nitro to amine reduction detected at depth {depth}")

                    # Check for azide reduction
                    elif any(checker.check_fg("Azide", r) for r in reactants_smiles):
                        if checker.check_reaction("Azide to amine reduction (Staudinger)", rsmi):
                            print(f"Azide to amine reduction detected at depth {depth}")

                    # Check for nitrile reduction
                    elif any(checker.check_fg("Nitrile", r) for r in reactants_smiles):
                        if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                            print(f"Nitrile to amine reduction detected at depth {depth}")

                    # Check for amide reduction
                    elif any(checker.check_fg("Primary amide", r) for r in reactants_smiles) or any(
                        checker.check_fg("Secondary amide", r) for r in reactants_smiles
                    ):
                        if checker.check_reaction(
                            "Reduction of primary amides to amines", rsmi
                        ) or checker.check_reaction(
                            "Reduction of secondary amides to amines", rsmi
                        ):
                            print(f"Amide to amine reduction detected at depth {depth}")
                    else:
                        print(f"Amine formation detected at depth {depth}")

                    # Add to transformations regardless of source
                    transformations.append(("amine", depth))

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort transformations by depth (retrosynthetic order)
    transformations.sort(key=lambda x: x[1])
    sequence = [t[0] for t in transformations]

    print(f"Detected functional group sequence: {sequence}")

    # Define target sequence in retrosynthetic order
    target_sequence = ["amine", "heterocycle", "aldehyde", "vinyl"]

    # Check if the sequence contains our target transformations in the correct order
    if len(sequence) < len(target_sequence):
        print("Sequence too short to match target")
        return False

    # Check if all elements of target_sequence appear in sequence in the correct order
    # using a subsequence matching approach
    i, j = 0, 0  # i for target_sequence, j for sequence
    while i < len(target_sequence) and j < len(sequence):
        if target_sequence[i] == sequence[j]:
            print(f"Found {target_sequence[i]} at position {j}")
            i += 1
        j += 1

    # If we found all elements of target_sequence
    if i == len(target_sequence):
        print("Sequence matches target in correct order")
        return True

    print("Sequence elements not found in correct order")
    return False
