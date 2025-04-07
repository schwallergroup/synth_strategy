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
    Detects if the route uses a late-stage reductive amination to couple fragments.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Examining reaction at depth {depth}: {rsmi}")

            # Check if this is a late-stage reaction (depth <= 2)
            if depth <= 2:
                # Check if this is a reductive amination reaction
                is_reductive_amination = (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                )

                print(
                    f"Is reductive amination by reaction check: {is_reductive_amination}"
                )

                # Extract reactants to verify fragment coupling
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Need at least 2 reactants for fragment coupling
                    if len(reactants) >= 2:
                        # Check if one reactant has an amine and another has a carbonyl
                        has_amine = False
                        has_carbonyl = False
                        amine_reactants = []
                        carbonyl_reactants = []

                        for reactant in reactants:
                            # Check for various amine types
                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Aniline", reactant)
                            ):
                                has_amine = True
                                amine_reactants.append(reactant)
                                print(f"Found amine reactant: {reactant}")

                            # Check for carbonyl or alcohol groups
                            if (
                                checker.check_fg("Aldehyde", reactant)
                                or checker.check_fg("Ketone", reactant)
                                or checker.check_fg("Primary alcohol", reactant)
                                or checker.check_fg("Secondary alcohol", reactant)
                                or checker.check_fg("Formaldehyde", reactant)
                            ):
                                has_carbonyl = True
                                carbonyl_reactants.append(reactant)
                                print(f"Found carbonyl/alcohol reactant: {reactant}")

                        # If we have both amine and carbonyl reactants, it's coupling fragments
                        if has_amine and has_carbonyl:
                            print(
                                f"Confirmed fragment coupling via reductive amination"
                            )
                            found_pattern = True

                        # If reaction check failed but we have the right reactants, try to detect manually
                        elif not is_reductive_amination and has_amine and has_carbonyl:
                            # Look for N-C bond formation pattern
                            print("Checking for N-C bond formation pattern manually")

                            # Check if product has a tertiary or secondary amine that wasn't in reactants
                            if checker.check_fg(
                                "Tertiary amine", product
                            ) or checker.check_fg("Secondary amine", product):

                                # If no tertiary amine in reactants but present in product, likely reductive amination
                                if not any(
                                    checker.check_fg("Tertiary amine", r)
                                    for r in reactants
                                ):
                                    print("Detected N-C bond formation pattern")
                                    found_pattern = True

                        # Special case: Check for specific pattern in the reaction at depth 1
                        if depth == 1 and not found_pattern:
                            # Check if this looks like a reductive amination based on reactants and products
                            if any(
                                checker.check_fg("Aldehyde", r) for r in reactants
                            ) and any(
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Secondary amine", r)
                                for r in reactants
                            ):

                                # Check if product has a new C-N bond
                                if checker.check_fg(
                                    "Tertiary amine", product
                                ) or checker.check_fg("Secondary amine", product):
                                    print(
                                        "Detected reductive amination pattern at depth 1"
                                    )
                                    found_pattern = True

                except Exception as e:
                    print(f"Error analyzing reactants: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return found_pattern
