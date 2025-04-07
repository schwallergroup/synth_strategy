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
    This function detects if the synthesis includes a late-stage oxidation,
    particularly focusing on alcohol to carboxylic acid or aldehyde transformations
    in the first half of the synthesis depth.
    """
    late_oxidation_found = False
    max_depth = 0
    oxidation_depths = []

    # First pass: determine maximum depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass: identify oxidation reactions
    def find_oxidations(node, depth=0):
        nonlocal oxidation_depths

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for oxidation reactions using the checker function
            oxidation_reactions = [
                "Oxidation of alcohol to carboxylic acid",
                "Oxidation of aldehydes to carboxylic acids",
                "Oxidation of alcohol to aldehyde",
                "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
            ]

            for reaction_type in oxidation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Oxidation reaction detected: {reaction_type} at depth {depth}")
                    oxidation_depths.append(depth)
                    break

            # If no specific reaction type matched, check for functional group transformation
            if not oxidation_depths or depth not in oxidation_depths:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alcohol in reactants and carbonyl in product
                alcohol_in_reactants = any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants
                )

                aldehyde_in_product = checker.check_fg("Aldehyde", product)
                carboxylic_acid_in_product = checker.check_fg("Carboxylic acid", product)

                if alcohol_in_reactants and (aldehyde_in_product or carboxylic_acid_in_product):
                    print(f"Functional group oxidation detected at depth {depth}")
                    oxidation_depths.append(depth)

        # Traverse children
        for child in node.get("children", []):
            find_oxidations(child, depth + 1)

    # Execute passes
    find_max_depth(route)
    find_oxidations(route)

    # Determine if any oxidations are in the first half (late stage)
    if oxidation_depths:
        half_depth = max_depth / 2
        for depth in oxidation_depths:
            if depth <= half_depth:  # First half of synthesis (late stage)
                late_oxidation_found = True
                print(f"Late-stage oxidation found at depth {depth} (max depth: {max_depth})")
                break

    return late_oxidation_found
