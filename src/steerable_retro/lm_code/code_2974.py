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
    Detects if the route employs a sequence of redox manipulations
    (reduction followed by oxidation or vice versa).
    """
    # Track redox reactions with their depths
    redox_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")
                product = product_part

                # In retrosynthesis, the product is what we start with and reactants are what we're making
                # So for redox reactions, we need to check if product→reactants is a reduction or oxidation

                # Check for reduction reactions (in forward direction)
                # In retrosynthesis, these are oxidations (reactants → product)
                if (
                    checker.check_reaction("Reduction of aldehydes and ketones to alcohols", rsmi)
                    or checker.check_reaction("Reduction of ester to primary alcohol", rsmi)
                    or checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                    or checker.check_reaction(
                        "Reduction of carboxylic acid to primary alcohol", rsmi
                    )
                    or checker.check_reaction("Reduction of nitrile to amine", rsmi)
                    or checker.check_reaction("Reduction of primary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of secondary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of tertiary amides to amines", rsmi)
                ):
                    print(f"Found oxidation reaction at depth {depth}: {rsmi}")
                    redox_reactions.append(("oxidation", depth))

                # Check for oxidation reactions (in forward direction)
                # In retrosynthesis, these are reductions (reactants → product)
                elif (
                    checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                    or checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                    )
                    or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of ketone to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of amide to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of alkene to aldehyde", rsmi)
                    or checker.check_reaction("Oxidation of alcohol and aldehyde to ester", rsmi)
                    or checker.check_reaction("Oxidation of boronic acids", rsmi)
                    or checker.check_reaction("Oxidation of boronic esters", rsmi)
                    or checker.check_reaction("Oxidative esterification of primary alcohols", rsmi)
                ):
                    print(f"Found reduction reaction at depth {depth}: {rsmi}")
                    redox_reactions.append(("reduction", depth))

                # If no specific reaction type is found, check for functional group changes
                else:
                    # Check for functional groups in product and reactants
                    product_has_aldehyde = checker.check_fg("Aldehyde", product)
                    product_has_ketone = checker.check_fg("Ketone", product)
                    product_has_ester = checker.check_fg("Ester", product)
                    product_has_carboxylic = checker.check_fg("Carboxylic acid", product)
                    product_has_nitrile = checker.check_fg("Nitrile", product)
                    product_has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )
                    product_has_alcohol = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    )
                    product_has_amine = (
                        checker.check_fg("Primary amine", product)
                        or checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                    )

                    reactants_have_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                    )
                    reactants_have_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        for r in reactants
                    )
                    reactants_have_aldehyde = any(
                        checker.check_fg("Aldehyde", r) for r in reactants
                    )
                    reactants_have_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
                    reactants_have_ester = any(checker.check_fg("Ester", r) for r in reactants)
                    reactants_have_carboxylic = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    reactants_have_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants)
                    reactants_have_amide = any(
                        checker.check_fg("Primary amide", r)
                        or checker.check_fg("Secondary amide", r)
                        or checker.check_fg("Tertiary amide", r)
                        for r in reactants
                    )

                    # In retrosynthesis:
                    # Higher oxidation state → lower oxidation state = reduction
                    # Lower oxidation state → higher oxidation state = oxidation

                    # Check for reduction in retrosynthesis:
                    # Product (higher oxidation) → reactants (lower oxidation)
                    if (
                        product_has_aldehyde
                        or product_has_ketone
                        or product_has_ester
                        or product_has_carboxylic
                        or product_has_nitrile
                        or product_has_amide
                    ) and (reactants_have_alcohol or reactants_have_amine):
                        print(f"Found reduction pattern at depth {depth}: {rsmi}")
                        redox_reactions.append(("reduction", depth))

                    # Check for oxidation in retrosynthesis:
                    # Product (lower oxidation) → reactants (higher oxidation)
                    elif (product_has_alcohol or product_has_amine) and (
                        reactants_have_aldehyde
                        or reactants_have_ketone
                        or reactants_have_ester
                        or reactants_have_carboxylic
                        or reactants_have_nitrile
                        or reactants_have_amide
                    ):
                        print(f"Found oxidation pattern at depth {depth}: {rsmi}")
                        redox_reactions.append(("oxidation", depth))

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if there's a sequence of redox manipulations
    has_redox_sequence = False

    # Sort by depth to make checking adjacent steps easier
    redox_reactions.sort(key=lambda x: x[1])
    print(f"Redox reactions found: {redox_reactions}")

    # Check for redox steps of different types with at most one step in between
    for i in range(len(redox_reactions) - 1):
        current_type, current_depth = redox_reactions[i]
        next_type, next_depth = redox_reactions[i + 1]

        if current_type != next_type and abs(current_depth - next_depth) <= 2:
            print(
                f"Found redox sequence: {current_type} at depth {current_depth}, {next_type} at depth {next_depth}"
            )
            has_redox_sequence = True
            break

    return has_redox_sequence
