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
    Detects a synthesis that includes a reductive amination step to form a tertiary amine.
    """
    has_reductive_amination = False

    def dfs_traverse(node):
        nonlocal has_reductive_amination

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Checking reaction: {rsmi}")

            # Check if this is a reductive amination reaction
            is_reductive_amination = (
                checker.check_reaction("Reductive amination with aldehyde", rsmi)
                or checker.check_reaction("Reductive amination with ketone", rsmi)
                or checker.check_reaction("Reductive amination with alcohol", rsmi)
            )

            print(f"Is reductive amination reaction: {is_reductive_amination}")

            # If not directly identified as reductive amination, check for characteristic patterns
            if not is_reductive_amination:
                # Check for aldehyde/ketone + amine → tertiary amine pattern
                has_aldehyde = any(
                    checker.check_fg("Aldehyde", r) for r in reactants_smiles
                )
                has_ketone = any(
                    checker.check_fg("Ketone", r) for r in reactants_smiles
                )
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants_smiles
                )
                has_secondary_amine = any(
                    checker.check_fg("Secondary amine", r) for r in reactants_smiles
                )
                has_tertiary_amine_product = checker.check_fg(
                    "Tertiary amine", product_smiles
                )

                print(f"Manual check - Aldehyde: {has_aldehyde}, Ketone: {has_ketone}")
                print(
                    f"Manual check - Primary amine: {has_primary_amine}, Secondary amine: {has_secondary_amine}"
                )
                print(
                    f"Manual check - Tertiary amine in product: {has_tertiary_amine_product}"
                )

                # Check for reducing agents in the reagents section
                reagents = rsmi.split(">")[1].split(".")
                reducing_agents = [
                    "NaBH4",
                    "NaBH3CN",
                    "NaCNBH3",
                    "BC#N",
                    "B",
                    "H",
                    "[HH]",
                    "H2",
                ]
                has_reducing_agent = any(agent in reagents for agent in reducing_agents)

                print(f"Has reducing agent: {has_reducing_agent}")

                # If we have carbonyl + amine → tertiary amine + reducing agent, it's likely reductive amination
                if (
                    (has_aldehyde or has_ketone)
                    and (has_primary_amine or has_secondary_amine)
                    and has_tertiary_amine_product
                    and has_reducing_agent
                ):
                    is_reductive_amination = True
                    print("Identified as reductive amination through pattern matching")

            if is_reductive_amination:
                print("Detected reductive amination reaction")

                # Verify reactants contain aldehyde/ketone and amine
                has_aldehyde = any(
                    checker.check_fg("Aldehyde", r) for r in reactants_smiles
                )
                has_ketone = any(
                    checker.check_fg("Ketone", r) for r in reactants_smiles
                )
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants_smiles
                )
                has_secondary_amine = any(
                    checker.check_fg("Secondary amine", r) for r in reactants_smiles
                )

                # Verify product contains tertiary amine
                has_tertiary_amine = checker.check_fg("Tertiary amine", product_smiles)

                print(f"Aldehyde: {has_aldehyde}, Ketone: {has_ketone}")
                print(
                    f"Primary amine: {has_primary_amine}, Secondary amine: {has_secondary_amine}"
                )
                print(f"Tertiary amine in product: {has_tertiary_amine}")

                # Confirm reductive amination to form tertiary amine
                if (
                    (has_aldehyde or has_ketone)
                    and (has_primary_amine or has_secondary_amine)
                    and has_tertiary_amine
                ):
                    has_reductive_amination = True
                    print("Confirmed reductive amination to form tertiary amine")
                # Special case: check if secondary amine forms tertiary amine
                elif has_secondary_amine and has_tertiary_amine:
                    has_reductive_amination = True
                    print(
                        "Confirmed reductive amination: secondary amine to tertiary amine"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_reductive_amination
