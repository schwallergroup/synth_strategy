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
    This function detects redox sequences in a synthetic route.
    It looks for reduction followed by oxidation or vice versa within a few steps.
    """
    # Store reactions by depth for sequence analysis
    reactions_by_depth = {}
    redox_sequence_found = False

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Store reaction at its depth
            reactions_by_depth[depth] = node["metadata"]["rsmi"]

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # First pass: collect reactions by depth
    dfs_traverse(route)

    # Define oxidation and reduction reaction types to check
    oxidation_reactions = [
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        "Oxidation of alcohol to carboxylic acid",
        "Oxidation of alkene to carboxylic acid",
        "Oxidation of ketone to carboxylic acid",
        "Oxidation of nitrile to carboxylic acid",
        "Oxidation of alcohol and aldehyde to ester",
        "Oxidation of boronic acids",
        "Oxidation of boronic esters",
        "Oxidative esterification of primary alcohols",
        "Oxidative Heck reaction",
        "Aromatic hydroxylation",
    ]

    reduction_reactions = [
        "Reduction of aldehydes and ketones to alcohols",
        "Reduction of ester to primary alcohol",
        "Reduction of ketone to secondary alcohol",
        "Reduction of carboxylic acid to primary alcohol",
        "Reduction of nitrile to amide",
        "Reduction of nitrile to amine",
        "Reduction of primary amides to amines",
        "Reduction of secondary amides to amines",
        "Reduction of tertiary amides to amines",
        "Hydrogenation (double to single)",
        "Hydrogenation (triple to double)",
        "Arene hydrogenation",
    ]

    # Check consecutive depths for redox patterns
    depths = sorted(reactions_by_depth.keys())

    for i in range(len(depths) - 1):
        # Check for consecutive or near-consecutive reactions (allow 1 step gap)
        for j in range(i + 1, min(i + 3, len(depths))):
            current_depth = depths[i]
            next_depth = depths[j]

            if next_depth - current_depth > 2:
                continue  # Skip if too far apart

            current_rxn = reactions_by_depth[current_depth]
            next_rxn = reactions_by_depth[next_depth]

            # Check for oxidation followed by reduction
            is_current_oxidation = any(
                checker.check_reaction(rxn, current_rxn) for rxn in oxidation_reactions
            )
            is_next_reduction = any(
                checker.check_reaction(rxn, next_rxn) for rxn in reduction_reactions
            )

            if is_current_oxidation and is_next_reduction:
                print(
                    f"Redox sequence detected: Oxidation at depth {current_depth} followed by Reduction at depth {next_depth}"
                )
                redox_sequence_found = True
                break

            # Check for reduction followed by oxidation
            is_current_reduction = any(
                checker.check_reaction(rxn, current_rxn) for rxn in reduction_reactions
            )
            is_next_oxidation = any(
                checker.check_reaction(rxn, next_rxn) for rxn in oxidation_reactions
            )

            if is_current_reduction and is_next_oxidation:
                print(
                    f"Redox sequence detected: Reduction at depth {current_depth} followed by Oxidation at depth {next_depth}"
                )
                redox_sequence_found = True
                break

    # If no reaction-based detection, try functional group-based detection
    if not redox_sequence_found:
        for i in range(len(depths) - 1):
            # Check for consecutive or near-consecutive reactions (allow 1 step gap)
            for j in range(i + 1, min(i + 3, len(depths))):
                current_depth = depths[i]
                next_depth = depths[j]

                if next_depth - current_depth > 2:
                    continue  # Skip if too far apart

                current_rxn = reactions_by_depth[current_depth]
                next_rxn = reactions_by_depth[next_depth]

                # Extract products and reactants
                current_reactants = current_rxn.split(">")[0].split(".")
                current_product = current_rxn.split(">")[-1]

                next_reactants = next_rxn.split(">")[0].split(".")
                next_product = next_rxn.split(">")[-1]

                # Check for alcohol → aldehyde/ketone/acid → alcohol pattern
                if checker.check_fg(
                    "Primary alcohol", current_product
                ) or checker.check_fg("Secondary alcohol", current_product):
                    # Check if next reaction oxidizes the alcohol
                    if (
                        checker.check_fg("Aldehyde", next_product)
                        or checker.check_fg("Ketone", next_product)
                        or checker.check_fg("Carboxylic acid", next_product)
                    ):

                        # Look ahead for reduction back to alcohol
                        if j + 1 < len(depths) and depths[j + 1] - next_depth <= 2:
                            next_next_rxn = reactions_by_depth[depths[j + 1]]
                            next_next_product = next_next_rxn.split(">")[-1]

                            if checker.check_fg(
                                "Primary alcohol", next_next_product
                            ) or checker.check_fg(
                                "Secondary alcohol", next_next_product
                            ):
                                print(
                                    f"Redox sequence detected: Alcohol → Carbonyl → Alcohol at depths {current_depth}, {next_depth}, {depths[j+1]}"
                                )
                                redox_sequence_found = True
                                break

                # Check for aldehyde → alcohol → carbonyl pattern
                if checker.check_fg("Aldehyde", current_product):
                    # Check if next reaction reduces the aldehyde
                    if checker.check_fg("Primary alcohol", next_product):
                        # Look ahead for oxidation back to carbonyl
                        if j + 1 < len(depths) and depths[j + 1] - next_depth <= 2:
                            next_next_rxn = reactions_by_depth[depths[j + 1]]
                            next_next_product = next_next_rxn.split(">")[-1]

                            if (
                                checker.check_fg("Aldehyde", next_next_product)
                                or checker.check_fg("Ketone", next_next_product)
                                or checker.check_fg(
                                    "Carboxylic acid", next_next_product
                                )
                                or checker.check_fg("Ester", next_next_product)
                            ):
                                print(
                                    f"Redox sequence detected: Aldehyde → Alcohol → Carbonyl at depths {current_depth}, {next_depth}, {depths[j+1]}"
                                )
                                redox_sequence_found = True
                                break

    return redox_sequence_found
