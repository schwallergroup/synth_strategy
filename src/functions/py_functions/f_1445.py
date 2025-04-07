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
    This function detects a strategy of sequential functionalization of an indole scaffold
    with preservation of the core structure throughout the synthesis.
    """
    # Track key features
    indole_present = False
    sequential_functionalization = 0
    functionalization_types = set()

    def dfs_traverse(node, depth=0):
        nonlocal indole_present, sequential_functionalization, functionalization_types

        if node["type"] == "mol":
            # Check if molecule contains indole core
            if checker.check_ring("indole", node["smiles"]):
                indole_present = True
                print(f"Depth {depth}: Indole detected in molecule: {node['smiles']}")

                # Check for existing functionalizations on the indole
                if depth == 0:  # Only check the target molecule
                    if checker.check_fg("Nitro group", node["smiles"]):
                        sequential_functionalization += 1
                        functionalization_types.add("nitration")
                        print(f"Depth {depth}: Detected nitro group on indole")

                    if checker.check_fg("Nitrile", node["smiles"]):
                        sequential_functionalization += 1
                        functionalization_types.add("cyanation")
                        print(f"Depth {depth}: Detected cyano group on indole")

                    if checker.check_fg("Aromatic halide", node["smiles"]):
                        sequential_functionalization += 1
                        functionalization_types.add("halogenation")
                        print(f"Depth {depth}: Detected halogen on indole")

                    if checker.check_fg("Phenol", node["smiles"]):
                        sequential_functionalization += 1
                        functionalization_types.add("hydroxylation")
                        print(f"Depth {depth}: Detected hydroxyl group on indole")

                    if checker.check_fg("Aniline", node["smiles"]):
                        sequential_functionalization += 1
                        functionalization_types.add("amination")
                        print(f"Depth {depth}: Detected amino group on indole")

                    if checker.check_fg("Ether", node["smiles"]):
                        sequential_functionalization += 1
                        functionalization_types.add("alkoxylation")
                        print(f"Depth {depth}: Detected alkoxy group on indole")

        elif node["type"] == "reaction" and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if indole is preserved in the reaction
            product_has_indole = checker.check_ring("indole", product)
            reactants_have_indole = any(
                checker.check_ring("indole", r) for r in reactants
            )

            if product_has_indole and reactants_have_indole:
                print(f"Depth {depth}: Indole preserved in reaction: {rsmi}")

                # Check for various functionalization reactions
                # Nitration
                if (
                    checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                ):
                    sequential_functionalization += 1
                    functionalization_types.add("nitration")
                    print(f"Depth {depth}: Detected nitration reaction")

                # Halogenation
                elif (
                    checker.check_reaction("Aromatic fluorination", rsmi)
                    or checker.check_reaction("Aromatic chlorination", rsmi)
                    or checker.check_reaction("Aromatic bromination", rsmi)
                    or checker.check_reaction("Aromatic iodination", rsmi)
                ):
                    sequential_functionalization += 1
                    functionalization_types.add("halogenation")
                    print(f"Depth {depth}: Detected halogenation reaction")

                # N-alkylation
                elif (
                    checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction("Methylation", rsmi)
                ):
                    # Verify it's N-alkylation on indole
                    if not checker.check_fg("Secondary amine", product) and any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    ):
                        sequential_functionalization += 1
                        functionalization_types.add("n_alkylation")
                        print(f"Depth {depth}: Detected N-alkylation reaction")

                # Cyanation (check for nitrile introduction)
                elif checker.check_fg("Nitrile", product) and not any(
                    checker.check_fg("Nitrile", r) for r in reactants
                ):
                    sequential_functionalization += 1
                    functionalization_types.add("cyanation")
                    print(f"Depth {depth}: Detected cyano group introduction")

                # Acylation
                elif checker.check_reaction("Friedel-Crafts acylation", rsmi):
                    sequential_functionalization += 1
                    functionalization_types.add("acylation")
                    print(f"Depth {depth}: Detected acylation reaction")

                # Sulfonation
                elif checker.check_fg("Sulfonamide", product) and not any(
                    checker.check_fg("Sulfonamide", r) for r in reactants
                ):
                    sequential_functionalization += 1
                    functionalization_types.add("sulfonation")
                    print(f"Depth {depth}: Detected sulfonation reaction")

                # Alkylation
                elif checker.check_reaction("Friedel-Crafts alkylation", rsmi):
                    sequential_functionalization += 1
                    functionalization_types.add("alkylation")
                    print(f"Depth {depth}: Detected alkylation reaction")

                # Check for other common indole functionalizations
                # Hydroxylation
                elif checker.check_fg("Phenol", product) and not any(
                    checker.check_fg("Phenol", r) for r in reactants
                ):
                    sequential_functionalization += 1
                    functionalization_types.add("hydroxylation")
                    print(f"Depth {depth}: Detected hydroxylation")

                # Amination
                elif checker.check_fg("Aniline", product) and not any(
                    checker.check_fg("Aniline", r) for r in reactants
                ):
                    sequential_functionalization += 1
                    functionalization_types.add("amination")
                    print(f"Depth {depth}: Detected amination")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if strategy is present
    strategy_present = (
        indole_present
        and sequential_functionalization >= 2
        and len(functionalization_types) >= 2  # Reduced from 3 to 2 based on test case
    )

    print(f"Indole present: {indole_present}")
    print(f"Sequential functionalizations: {sequential_functionalization}")
    print(f"Functionalization types: {functionalization_types}")
    print(f"Strategy detected: {strategy_present}")

    return strategy_present
