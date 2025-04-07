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
    This function detects if the synthesis uses late-stage reductive amination
    to form a tertiary amine in the final steps.
    """
    # Initialize tracking variables
    has_reductive_amination = False
    reductive_amination_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal has_reductive_amination, reductive_amination_depth

        # Adjust depth calculation - only count reactions for depth
        current_depth = depth

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            try:
                # Check if this is a reductive amination reaction
                is_reductive_amination = (
                    checker.check_reaction("reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("reductive amination with ketone", rsmi)
                    or checker.check_reaction("reductive amination with alcohol", rsmi)
                    or checker.check_reaction(
                        "Eschweiler-Clarke Primary Amine Methylation", rsmi
                    )
                    or checker.check_reaction(
                        "Eschweiler-Clarke Secondary Amine Methylation", rsmi
                    )
                    or checker.check_reaction(
                        "reductive methylation of primary amine with formaldehyde", rsmi
                    )
                )

                print(
                    f"Reaction check results: reductive_amination_aldehyde={checker.check_reaction('reductive amination with aldehyde', rsmi)}, "
                    f"reductive_amination_ketone={checker.check_reaction('reductive amination with ketone', rsmi)}, "
                    f"reductive_amination_alcohol={checker.check_reaction('reductive amination with alcohol', rsmi)}, "
                    f"eschweiler_clarke_primary={checker.check_reaction('Eschweiler-Clarke Primary Amine Methylation', rsmi)}, "
                    f"eschweiler_clarke_secondary={checker.check_reaction('Eschweiler-Clarke Secondary Amine Methylation', rsmi)}, "
                    f"reductive_methylation={checker.check_reaction('reductive methylation of primary amine with formaldehyde', rsmi)}"
                )

                if is_reductive_amination:
                    print(f"Found reductive amination reaction at depth {depth}")

                # Check if tertiary amine is formed in the product but not present in reactants
                has_tertiary_amine_product = checker.check_fg("Tertiary amine", product)
                tertiary_amine_in_reactants = any(
                    checker.check_fg("Tertiary amine", r) for r in reactants
                )

                if has_tertiary_amine_product:
                    print(f"Product has tertiary amine")
                if tertiary_amine_in_reactants:
                    print(f"Reactants already have tertiary amine")

                # Verify we have either aldehyde or ketone in reactants
                has_aldehyde = any(
                    checker.check_fg("Aldehyde", r)
                    or checker.check_fg("Formaldehyde", r)
                    for r in reactants
                )
                has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
                has_alcohol = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants
                )

                # Verify we have primary or secondary amine in reactants
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants
                )
                has_secondary_amine = any(
                    checker.check_fg("Secondary amine", r) for r in reactants
                )

                print(
                    f"Reactants have: aldehyde={has_aldehyde}, ketone={has_ketone}, alcohol={has_alcohol}, "
                    f"primary_amine={has_primary_amine}, secondary_amine={has_secondary_amine}"
                )

                # Check for formaldehyde specifically
                has_formaldehyde = any(
                    checker.check_fg("Formaldehyde", r) for r in reactants
                )
                if has_formaldehyde:
                    print(f"Formaldehyde detected in reactants")

                # Modified condition to better detect reductive amination
                if (
                    (
                        is_reductive_amination
                        or (
                            has_formaldehyde
                            and has_primary_amine
                            and has_tertiary_amine_product
                        )
                    )
                    and has_tertiary_amine_product
                    and not tertiary_amine_in_reactants
                    and (has_aldehyde or has_ketone or has_alcohol or has_formaldehyde)
                    and (has_primary_amine or has_secondary_amine)
                ):
                    print(
                        f"Confirmed reductive amination forming tertiary amine at depth {depth}"
                    )
                    has_reductive_amination = True
                    reductive_amination_depth = min(reductive_amination_depth, depth)
            except Exception as e:
                print(f"Error processing reaction SMILES: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if reductive amination is detected in the late stage (depth ≤ 2)
    # In retrosynthesis, depth 0 is target, depth 1 is first reaction, depth 2 is reactants of first reaction
    late_stage = reductive_amination_depth <= 2
    print(f"Reductive amination depth: {reductive_amination_depth}")
    print(
        f"Late-stage reductive amination detected: {has_reductive_amination and late_stage}"
    )
    return has_reductive_amination and late_stage
