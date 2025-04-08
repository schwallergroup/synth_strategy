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
    This function detects if the synthesis uses a late-stage reductive amination strategy
    (conversion of a ketone or aldehyde to an amine in the final steps).
    """
    found_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal found_reductive_amination

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is a reductive amination reaction
                is_reductive_amination = (
                    checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                    or checker.check_reaction("reductive amination", rsmi)
                )

                if is_reductive_amination:
                    print(f"Detected reductive amination reaction at depth {depth}")

                    # Extract reactants and product to verify functional group transformation
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check for ketone, aldehyde, or alcohol in reactants
                    has_carbonyl = any(
                        checker.check_fg("Ketone", r)
                        or checker.check_fg("Aldehyde", r)
                        or checker.check_fg(
                            "Alcohol", r
                        )  # Include alcohols which can be oxidized in situ
                        for r in reactants
                        if r
                    )

                    # Check for ammonia in reactants
                    has_ammonia = any("NH3" in r for r in reactants)

                    # Check for amine in product
                    has_amine = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    )

                    print(f"Carbonyl/alcohol in reactants: {has_carbonyl}")
                    print(f"Ammonia in reactants: {has_ammonia}")
                    print(f"Amine in product: {has_amine}")

                    if (has_carbonyl and has_amine) or (has_carbonyl and has_ammonia):
                        print(f"Found late-stage reductive amination at depth {depth}")
                        found_reductive_amination = True
                else:
                    # Try to detect reductive amination by examining functional group changes
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for ketone, aldehyde, or alcohol in reactants
                    has_carbonyl = any(
                        checker.check_fg("Ketone", r)
                        or checker.check_fg("Aldehyde", r)
                        or checker.check_fg("Alcohol", r)
                        for r in reactants
                        if r
                    )

                    # Check for nitrogen source in reactants
                    has_nitrogen_source = any(
                        "NH3" in r
                        or checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                        if r
                    )

                    # Check for amine in product
                    has_amine = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    )

                    # Check for reducing agents or conditions
                    has_reducing_conditions = any(
                        "B" in r
                        or "NaBH" in r
                        or "LiAlH" in r
                        or "Ti" in r
                        or "H2" in r
                        or "[H]" in r
                        for r in reactants
                        if r
                    )

                    if (
                        has_carbonyl
                        and has_nitrogen_source
                        and has_amine
                        and has_reducing_conditions
                    ):
                        print(
                            f"Found late-stage reductive amination by functional group analysis at depth {depth}"
                        )
                        found_reductive_amination = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: {found_reductive_amination}")
    return found_reductive_amination
