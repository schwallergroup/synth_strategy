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
    This function detects multiple redox transformations in a synthesis route,
    specifically looking for at least one reduction and one oxidation.
    """
    reductions = 0
    oxidations = 0

    # List of reduction reaction types
    reduction_reactions = [
        "Reduction of aldehydes and ketones to alcohols",
        "Reduction of carboxylic acid to primary alcohol",
        "Reduction of ester to primary alcohol",
        "Reduction of nitrile to amine",
        "Reduction of nitro groups to amines",
        "Reduction of primary amides to amines",
        "Reduction of secondary amides to amines",
        "Reduction of tertiary amides to amines",
        "Hydrogenation (double to single)",
        "Hydrogenation (triple to double)",
        "Arene hydrogenation",
        "Grignard from aldehyde to alcohol",
        "Grignard from ketone to alcohol",
    ]

    # List of oxidation reaction types
    oxidation_reactions = [
        "Oxidation of aldehydes to carboxylic acids",
        "Oxidation of alkene to carboxylic acid",
        "Oxidation of alcohol to carboxylic acid",
        "Oxidation of ketone to carboxylic acid",
        "Oxidation of nitrile to carboxylic acid",
        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
        "Oxidation of alkene to aldehyde",
        "Oxidative esterification of primary alcohols",
        "Oxidation of alcohol and aldehyde to ester",
        "Quinone formation",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal reductions, oxidations

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # In retrosynthesis, the product is the starting material and reactants are the targets
            # So we need to check if the reaction in forward direction is a redox transformation

            # Check for reduction reactions
            for reaction_type in reduction_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    reductions += 1
                    print(f"Reduction detected: {reaction_type}, total: {reductions}")
                    break

            # Check for oxidation reactions
            for reaction_type in oxidation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    oxidations += 1
                    print(f"Oxidation detected: {reaction_type}, total: {oxidations}")
                    break

            # If no specific reaction type was found, try to detect redox transformations by functional group changes
            if reductions == 0 or oxidations == 0:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for reduction patterns (in forward direction)
                if any(
                    checker.check_fg("Aldehyde", r) for r in reactants
                ) and checker.check_fg("Primary alcohol", product):
                    reductions += 1
                    print(
                        f"Reduction detected (Aldehyde → Primary alcohol), total: {reductions}"
                    )

                elif any(
                    checker.check_fg("Ketone", r) for r in reactants
                ) and checker.check_fg("Secondary alcohol", product):
                    reductions += 1
                    print(
                        f"Reduction detected (Ketone → Secondary alcohol), total: {reductions}"
                    )

                elif any(
                    checker.check_fg("Nitrile", r) for r in reactants
                ) and checker.check_fg("Primary amine", product):
                    reductions += 1
                    print(
                        f"Reduction detected (Nitrile → Primary amine), total: {reductions}"
                    )

                elif any(
                    checker.check_fg("Nitro group", r) for r in reactants
                ) and checker.check_fg("Primary amine", product):
                    reductions += 1
                    print(
                        f"Reduction detected (Nitro group → Primary amine), total: {reductions}"
                    )

                # Additional check for Grignard-type reductions
                elif any(
                    checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                    for r in reactants
                ) and (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                ):
                    reductions += 1
                    print(
                        f"Reduction detected (Carbonyl → Alcohol), total: {reductions}"
                    )

                # Check for oxidation patterns (in forward direction)
                if any(
                    checker.check_fg("Primary alcohol", r) for r in reactants
                ) and checker.check_fg("Aldehyde", product):
                    oxidations += 1
                    print(
                        f"Oxidation detected (Primary alcohol → Aldehyde), total: {oxidations}"
                    )

                elif any(
                    checker.check_fg("Secondary alcohol", r) for r in reactants
                ) and checker.check_fg("Ketone", product):
                    oxidations += 1
                    print(
                        f"Oxidation detected (Secondary alcohol → Ketone), total: {oxidations}"
                    )

                elif any(
                    checker.check_fg("Primary alcohol", r) for r in reactants
                ) and checker.check_fg("Carboxylic acid", product):
                    oxidations += 1
                    print(
                        f"Oxidation detected (Primary alcohol → Carboxylic acid), total: {oxidations}"
                    )

                elif any(
                    checker.check_fg("Aldehyde", r) for r in reactants
                ) and checker.check_fg("Carboxylic acid", product):
                    oxidations += 1
                    print(
                        f"Oxidation detected (Aldehyde → Carboxylic acid), total: {oxidations}"
                    )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    print("Starting traversal to find redox transformations...")
    dfs_traverse(route)

    print(f"Final count - Reductions: {reductions}, Oxidations: {oxidations}")
    # Return True if at least one reduction and one oxidation were found
    return reductions >= 1 and oxidations >= 1
