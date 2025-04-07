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
    This function detects if the synthetic route involves a late-stage reduction of a carbonyl group.
    """
    reduction_at_low_depth = False

    def dfs_traverse(node, depth=0):
        nonlocal reduction_at_low_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for carbonyl reduction in either direction

                # Forward direction: alcohol → carbonyl (oxidation)
                has_alcohol_reactant = any(
                    r
                    and (
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                    )
                    for r in reactants
                )

                has_carbonyl_product = product and (
                    checker.check_fg("Aldehyde", product)
                    or checker.check_fg("Ketone", product)
                    or checker.check_fg("Carboxylic acid", product)
                )

                # Reverse direction: carbonyl → alcohol (reduction)
                has_carbonyl_reactant = any(
                    r
                    and (
                        checker.check_fg("Aldehyde", r)
                        or checker.check_fg("Ketone", r)
                        or checker.check_fg("Carboxylic acid", r)
                    )
                    for r in reactants
                )

                has_alcohol_product = product and (
                    checker.check_fg("Primary alcohol", product)
                    or checker.check_fg("Secondary alcohol", product)
                    or checker.check_fg("Tertiary alcohol", product)
                    or checker.check_fg("Aromatic alcohol", product)
                )

                # Check for various oxidation/reduction reaction types
                is_oxidation = (
                    checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Oxidation of alcohol to carboxylic acid", rsmi
                    )
                    or checker.check_reaction("Oxidation of primary alcohols", rsmi)
                )

                is_reduction = (
                    checker.check_reaction(
                        "Reduction of aldehydes and ketones to alcohols", rsmi
                    )
                    or checker.check_reaction(
                        "Reduction of carboxylic acid to primary alcohol", rsmi
                    )
                    or checker.check_reaction(
                        "Reduction of ester to primary alcohol", rsmi
                    )
                )

                # Check for pattern matching if reaction type check fails
                pattern_match_oxidation = has_alcohol_reactant and has_carbonyl_product
                pattern_match_reduction = has_carbonyl_reactant and has_alcohol_product

                is_reduction_or_oxidation = (
                    is_oxidation
                    or is_reduction
                    or pattern_match_oxidation
                    or pattern_match_reduction
                )

                print(
                    f"Has alcohol in reactants: {has_alcohol_reactant}, Has carbonyl in product: {has_carbonyl_product}"
                )
                print(
                    f"Has carbonyl in reactants: {has_carbonyl_reactant}, Has alcohol in product: {has_alcohol_product}"
                )
                print(f"Is oxidation: {is_oxidation}, Is reduction: {is_reduction}")
                print(f"Is reduction or oxidation: {is_reduction_or_oxidation}")

                # Consider it a late-stage reduction if it occurs at a low depth (≤ 2)
                if is_reduction_or_oxidation and depth <= 2:
                    print(f"Late-stage carbonyl reduction detected at depth {depth}")
                    reduction_at_low_depth = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return reduction_at_low_depth
