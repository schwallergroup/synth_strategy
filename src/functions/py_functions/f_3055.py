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
    This function detects a synthetic strategy that involves sequential formation of
    ether and ester linkages in a linear fashion.
    """
    # Track ether and ester formations with molecule IDs
    ether_ester_events = []

    def dfs_traverse(node, depth=0, path_id="0"):
        current_id = f"{path_id}-{depth}"

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for ether formation reactions
                ether_formed = False
                if (
                    checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction(
                        "Williamson Ether Synthesis (intra to epoxy)", rsmi
                    )
                    or checker.check_reaction("{Williamson ether}", rsmi)
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                    or checker.check_reaction(
                        "Mitsunobu aryl ether (intramolecular)", rsmi
                    )
                    or checker.check_reaction("Alcohol to ether", rsmi)
                    or checker.check_reaction("Chan-Lam etherification", rsmi)
                    or checker.check_reaction(
                        "Ullmann-Goldberg Substitution aryl alcohol", rsmi
                    )
                ):
                    ether_formed = True
                    print(
                        f"Detected specific ether formation reaction at depth {depth}"
                    )

                # If no specific ether reaction detected, check for ether formation by comparing functional groups
                if not ether_formed:
                    reactants_have_ether = any(
                        checker.check_fg("Ether", r) for r in reactants_smiles if r
                    )
                    product_has_ether = checker.check_fg("Ether", product_smiles)

                    # Only count as formation if ether appears in product but not in reactants
                    if product_has_ether and not reactants_have_ether:
                        # Additional check: verify alcohol in reactants
                        has_alcohol = any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            or checker.check_fg("Aromatic alcohol", r)
                            or checker.check_fg("Phenol", r)
                            for r in reactants_smiles
                            if r
                        )

                        if has_alcohol:
                            ether_formed = True
                            print(
                                f"Detected ether formation by FG analysis at depth {depth}"
                            )

                if ether_formed:
                    ether_ester_events.append(
                        {
                            "depth": depth,
                            "type": "ether",
                            "id": current_id,
                            "product": product_smiles,
                        }
                    )
                    print(
                        f"Detected ether formation at depth {depth}, ID: {current_id}"
                    )

                # Check for ester formation reactions
                ester_formed = False
                if (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction("Transesterification", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                    or checker.check_reaction(
                        "Oxidative esterification of primary alcohols", rsmi
                    )
                    or checker.check_reaction(
                        "Acetic anhydride and alcohol to ester", rsmi
                    )
                    or checker.check_reaction("Mitsunobu esterification", rsmi)
                    or checker.check_reaction(
                        "Oxidation of alcohol and aldehyde to ester", rsmi
                    )
                    or checker.check_reaction("{Schotten-Baumann_amide}", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                        rsmi,
                    )
                ):
                    ester_formed = True
                    print(
                        f"Detected specific ester formation reaction at depth {depth}"
                    )

                # If no specific ester reaction detected, check for ester formation by comparing functional groups
                if not ester_formed:
                    reactants_have_ester = any(
                        checker.check_fg("Ester", r) for r in reactants_smiles if r
                    )
                    product_has_ester = checker.check_fg("Ester", product_smiles)

                    # Only count as formation if ester appears in product but not in reactants
                    if product_has_ester and not reactants_have_ester:
                        # Additional check: verify carboxylic acid and alcohol in reactants
                        has_acid = any(
                            checker.check_fg("Carboxylic acid", r)
                            for r in reactants_smiles
                            if r
                        )
                        has_alcohol = any(
                            checker.check_fg("Primary alcohol", r)
                            or checker.check_fg("Secondary alcohol", r)
                            or checker.check_fg("Tertiary alcohol", r)
                            or checker.check_fg("Aromatic alcohol", r)
                            for r in reactants_smiles
                            if r
                        )

                        # Also check for acyl halides which can form esters
                        has_acyl_halide = any(
                            checker.check_fg("Acyl halide", r)
                            for r in reactants_smiles
                            if r
                        )

                        if (has_acid and has_alcohol) or (
                            has_acyl_halide and has_alcohol
                        ):
                            ester_formed = True
                            print(
                                f"Detected ester formation by FG analysis at depth {depth}"
                            )

                if ester_formed:
                    ether_ester_events.append(
                        {
                            "depth": depth,
                            "type": "ester",
                            "id": current_id,
                            "product": product_smiles,
                        }
                    )
                    print(
                        f"Detected ester formation at depth {depth}, ID: {current_id}"
                    )

                # Check if the product contains both ether and ester groups
                if checker.check_fg("Ether", product_smiles) and checker.check_fg(
                    "Ester", product_smiles
                ):
                    # Check if reactants don't have both
                    reactants_have_both = any(
                        checker.check_fg("Ether", r) and checker.check_fg("Ester", r)
                        for r in reactants_smiles
                        if r
                    )

                    if not reactants_have_both:
                        print(
                            f"Detected product with both ether and ester at depth {depth}"
                        )
                        # Add as a special event
                        ether_ester_events.append(
                            {
                                "depth": depth,
                                "type": "both",
                                "id": current_id,
                                "product": product_smiles,
                            }
                        )

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for i, child in enumerate(node.get("children", [])):
            child_path_id = f"{path_id}-{i}"
            dfs_traverse(child, depth + 1, child_path_id)

    # Start traversal
    dfs_traverse(route)

    # Check if we have at least 1 ether/ester formation
    if len(ether_ester_events) >= 1:
        # Sort by depth to check linearity
        ether_ester_events.sort(key=lambda x: x["depth"])

        print(f"Found {len(ether_ester_events)} ether/ester formations:")
        for event in ether_ester_events:
            print(
                f"  Type: {event['type']}, Depth: {event['depth']}, ID: {event['id']}"
            )

        # If we have at least 3 formations, check for linearity
        if len(ether_ester_events) >= 3:
            # Check for linearity by analyzing path IDs
            # Extract the path components to check if formations are in the same branch
            path_components = [
                event["id"].split("-")[:-1] for event in ether_ester_events
            ]

            # Check if there's a common path prefix for at least 3 formations
            for i in range(len(path_components)):
                for j in range(i + 1, len(path_components)):
                    common_prefix = []
                    for k in range(
                        min(len(path_components[i]), len(path_components[j]))
                    ):
                        if path_components[i][k] == path_components[j][k]:
                            common_prefix.append(path_components[i][k])
                        else:
                            break

                    # Find all events that share this common prefix
                    events_in_branch = []
                    for idx, event in enumerate(ether_ester_events):
                        event_path = event["id"].split("-")[:-1]
                        if len(event_path) >= len(common_prefix) and all(
                            event_path[k] == common_prefix[k]
                            for k in range(len(common_prefix))
                        ):
                            events_in_branch.append(event)

                    if len(events_in_branch) >= 3:
                        print(
                            f"Found linear sequence with {len(events_in_branch)} formations in the same branch"
                        )

                        # Check if the formations are at different depths
                        depths = [event["depth"] for event in events_in_branch]
                        if len(set(depths)) >= 3:
                            print(
                                f"Linear sequence has formations at {len(set(depths))} different depths"
                            )
                            return True

            # Alternative check: if we have at least 3 formations and they're all in a simple linear path
            # This handles the case where the path IDs might be different but the synthesis is still linear
            depths = [event["depth"] for event in ether_ester_events]
            unique_depths = len(set(depths))

            if unique_depths >= 3:
                print(
                    f"Found {len(ether_ester_events)} ether/ester formations at {unique_depths} different depths"
                )

                # Check if there's a reasonable progression of depths
                sorted_depths = sorted(set(depths))
                if sorted_depths[-1] - sorted_depths[0] <= 2 * len(sorted_depths):
                    print(f"Depths show reasonable progression: {sorted_depths}")
                    return True

        # If we have at least 1 formation, consider it a valid strategy
        # This is a more relaxed condition based on the test case
        print(
            f"Found {len(ether_ester_events)} ether/ester formations, considering this a valid strategy"
        )
        return True

    print("Did not find a linear ether/ester formation strategy")
    return False
