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
    This function detects late-stage reductive amination for C-N bond formation.
    """
    reductive_amination_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_detected

        if (
            node["type"] == "reaction" and depth <= 2
        ):  # Check late-stage reactions (depth 0, 1, or 2)
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is a reductive amination reaction
                is_reductive_amination = (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                    or checker.check_reaction("Mignonac reaction", rsmi)
                    or checker.check_reaction("reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("reductive amination with ketone", rsmi)
                    or checker.check_reaction("reductive amination with alcohol", rsmi)
                )

                print(
                    f"Is reductive amination according to checker: {is_reductive_amination}"
                )

                # Extract reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for carbonyl compounds or alcohol in reactants
                aldehyde_found = any(checker.check_fg("Aldehyde", r) for r in reactants)
                ketone_found = any(checker.check_fg("Ketone", r) for r in reactants)
                formaldehyde_found = any(
                    checker.check_fg("Formaldehyde", r) for r in reactants
                )
                alcohol_found = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    or checker.check_fg("Aromatic alcohol", r)
                    for r in reactants
                )

                # Check for amine in reactants
                primary_amine_found = any(
                    checker.check_fg("Primary amine", r) for r in reactants
                )
                secondary_amine_found = any(
                    checker.check_fg("Secondary amine", r) for r in reactants
                )

                # Check for secondary/tertiary amine in product
                secondary_amine_in_product = checker.check_fg(
                    "Secondary amine", product
                )
                tertiary_amine_in_product = checker.check_fg("Tertiary amine", product)

                print(
                    f"Aldehyde: {aldehyde_found}, Ketone: {ketone_found}, Formaldehyde: {formaldehyde_found}, Alcohol: {alcohol_found}"
                )
                print(
                    f"Primary amine: {primary_amine_found}, Secondary amine: {secondary_amine_found}"
                )
                print(
                    f"Secondary amine in product: {secondary_amine_in_product}, Tertiary amine in product: {tertiary_amine_in_product}"
                )

                # Manual pattern check for reductive amination
                manual_check = False
                if (
                    (
                        aldehyde_found
                        or ketone_found
                        or formaldehyde_found
                        or alcohol_found
                    )
                    and (primary_amine_found or secondary_amine_found)
                    and (secondary_amine_in_product or tertiary_amine_in_product)
                ):
                    manual_check = True
                    print("Manual pattern check for reductive amination: PASSED")

                # Verify reductive amination
                if is_reductive_amination or manual_check:
                    print(f"Late-stage reductive amination confirmed at depth {depth}")
                    reductive_amination_detected = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return reductive_amination_detected
