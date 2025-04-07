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
    Detects an oxidation-reduction-oxidation pattern in the synthetic route
    """
    # Track oxidation and reduction steps
    redox_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthetic traversal, the product is the starting point
                # and reactants are what we're moving towards
                # For forward synthesis direction, we need to consider:
                # - Oxidation in retrosynthesis appears as reduction in forward synthesis
                # - Reduction in retrosynthesis appears as oxidation in forward synthesis

                # Check for oxidation reactions (in forward direction)
                if (
                    checker.check_reaction(
                        "Oxidation of aldehydes to carboxylic acids", rsmi
                    )
                    or checker.check_reaction(
                        "Oxidation of ketone to carboxylic acid", rsmi
                    )
                    or checker.check_reaction(
                        "Oxidation of alcohol to carboxylic acid", rsmi
                    )
                    or checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Oxidation of alkene to carboxylic acid", rsmi
                    )
                    or checker.check_reaction("Oxidation of alkene to aldehyde", rsmi)
                    or checker.check_reaction(
                        "Oxidation of nitrile to carboxylic acid", rsmi
                    )
                    or checker.check_reaction("Aromatic hydroxylation", rsmi)
                    or checker.check_reaction("Quinone formation", rsmi)
                    or checker.check_reaction("Oxidation of boronic acids", rsmi)
                    or checker.check_reaction("Oxidation of boronic esters", rsmi)
                    or checker.check_reaction(
                        "Oxidative esterification of primary alcohols", rsmi
                    )
                    or checker.check_reaction(
                        "Oxidation of alcohol and aldehyde to ester", rsmi
                    )
                ):
                    print(f"Oxidation reaction found at depth {depth}: {rsmi}")
                    redox_sequence.append(("oxidation", depth))

                # Check for reduction reactions (in forward direction)
                elif (
                    checker.check_reaction(
                        "Reduction of aldehydes and ketones to alcohols", rsmi
                    )
                    or checker.check_reaction(
                        "Reduction of ester to primary alcohol", rsmi
                    )
                    or checker.check_reaction(
                        "Reduction of ketone to secondary alcohol", rsmi
                    )
                    or checker.check_reaction(
                        "Reduction of carboxylic acid to primary alcohol", rsmi
                    )
                    or checker.check_reaction("Reduction of nitrile to amine", rsmi)
                    or checker.check_reaction(
                        "Reduction of primary amides to amines", rsmi
                    )
                    or checker.check_reaction(
                        "Reduction of secondary amides to amines", rsmi
                    )
                    or checker.check_reaction(
                        "Reduction of tertiary amides to amines", rsmi
                    )
                    or checker.check_reaction("Hydrogenation (double to single)", rsmi)
                    or checker.check_reaction("Hydrogenation (triple to double)", rsmi)
                    or checker.check_reaction("Arene hydrogenation", rsmi)
                    or checker.check_reaction(
                        "Reduction of nitro groups to amines", rsmi
                    )
                    or checker.check_reaction(
                        "Azide to amine reduction (Staudinger)", rsmi
                    )
                ):
                    print(f"Reduction reaction found at depth {depth}: {rsmi}")
                    redox_sequence.append(("reduction", depth))

                # If no specific reaction type matches, check for functional group changes
                # that indicate oxidation or reduction
                else:
                    # Check for oxidation by functional group changes (in forward direction)
                    # Reactants -> Product (forward synthesis)
                    if (
                        (
                            any(
                                checker.check_fg("Primary alcohol", r)
                                for r in reactants
                            )
                            and checker.check_fg("Aldehyde", product)
                        )
                        or (
                            any(
                                checker.check_fg("Secondary alcohol", r)
                                for r in reactants
                            )
                            and checker.check_fg("Ketone", product)
                        )
                        or (
                            any(checker.check_fg("Aldehyde", r) for r in reactants)
                            and checker.check_fg("Carboxylic acid", product)
                        )
                        or (
                            any(
                                checker.check_fg("Primary alcohol", r)
                                for r in reactants
                            )
                            and checker.check_fg("Carboxylic acid", product)
                        )
                        or (
                            any(
                                checker.check_fg("Secondary alcohol", r)
                                for r in reactants
                            )
                            and checker.check_fg("Carboxylic acid", product)
                        )
                        or (
                            any(checker.check_fg("Alkene", r) for r in reactants)
                            and (
                                checker.check_fg("Aldehyde", product)
                                or checker.check_fg("Ketone", product)
                                or checker.check_fg("Carboxylic acid", product)
                            )
                        )
                        or (
                            any(checker.check_fg("Alkyne", r) for r in reactants)
                            and checker.check_fg("Carboxylic acid", product)
                        )
                    ):
                        print(f"Oxidation by FG change found at depth {depth}: {rsmi}")
                        redox_sequence.append(("oxidation", depth))

                    # Check for reduction by functional group changes (in forward direction)
                    # Reactants -> Product (forward synthesis)
                    elif (
                        (
                            any(checker.check_fg("Aldehyde", r) for r in reactants)
                            and checker.check_fg("Primary alcohol", product)
                        )
                        or (
                            any(checker.check_fg("Ketone", r) for r in reactants)
                            and checker.check_fg("Secondary alcohol", product)
                        )
                        or (
                            any(
                                checker.check_fg("Carboxylic acid", r)
                                for r in reactants
                            )
                            and (
                                checker.check_fg("Aldehyde", product)
                                or checker.check_fg("Primary alcohol", product)
                            )
                        )
                        or (
                            any(checker.check_fg("Ester", r) for r in reactants)
                            and checker.check_fg("Primary alcohol", product)
                        )
                        or (
                            any(checker.check_fg("Alkyne", r) for r in reactants)
                            and checker.check_fg("Alkene", product)
                        )
                        or (
                            any(checker.check_fg("Alkene", r) for r in reactants)
                            and checker.check_fg("Alkane", product)
                        )
                        or (
                            any(checker.check_fg("Nitro group", r) for r in reactants)
                            and checker.check_fg("Primary amine", product)
                        )
                        or (
                            any(checker.check_fg("Nitrile", r) for r in reactants)
                            and checker.check_fg("Primary amine", product)
                        )
                        or (
                            any(checker.check_fg("Azide", r) for r in reactants)
                            and checker.check_fg("Primary amine", product)
                        )
                    ):
                        print(f"Reduction by FG change found at depth {depth}: {rsmi}")
                        redox_sequence.append(("reduction", depth))
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth to get chronological order in forward synthesis
    # Lower depth = later stage in synthesis (closer to final product)
    redox_sequence.sort(key=lambda x: x[1], reverse=True)

    # Extract just the reaction types
    redox_types = [r[0] for r in redox_sequence]
    print(f"Redox sequence in forward synthesis order: {redox_types}")

    # Check for oxidation-reduction-oxidation pattern
    pattern_found = False

    # Check for consecutive oxidation-reduction-oxidation pattern
    for i in range(len(redox_types) - 2):
        if (
            redox_types[i] == "oxidation"
            and redox_types[i + 1] == "reduction"
            and redox_types[i + 2] == "oxidation"
        ):
            pattern_found = True
            print(
                f"Found consecutive oxidation-reduction-oxidation pattern at positions {i}, {i+1}, {i+2}"
            )
            break

    # If not found, check for reduction-oxidation pattern (partial pattern)
    if not pattern_found and len(redox_types) >= 2:
        for i in range(len(redox_types) - 1):
            if redox_types[i] == "reduction" and redox_types[i + 1] == "oxidation":
                pattern_found = True
                print(
                    f"Found partial reduction-oxidation pattern at positions {i}, {i+1}"
                )
                break

    # If still not found, check for non-consecutive pattern
    if not pattern_found and len(redox_types) >= 3:
        # Find indices of each type
        oxidation_indices = [i for i, x in enumerate(redox_types) if x == "oxidation"]
        reduction_indices = [i for i, x in enumerate(redox_types) if x == "reduction"]

        # Check if we have at least 2 oxidations and 1 reduction in the right order
        if len(oxidation_indices) >= 2 and len(reduction_indices) >= 1:
            for red_idx in reduction_indices:
                # Find oxidations before and after this reduction
                before_ox = [i for i in oxidation_indices if i < red_idx]
                after_ox = [i for i in oxidation_indices if i > red_idx]

                if before_ox and after_ox:
                    pattern_found = True
                    print(
                        f"Found non-consecutive oxidation-reduction-oxidation pattern: oxidation at {before_ox[-1]}, reduction at {red_idx}, oxidation at {after_ox[0]}"
                    )
                    break

    print(f"Oxidation-reduction-oxidation pattern detected: {pattern_found}")

    return pattern_found
