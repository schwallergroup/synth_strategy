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
    Detects a strategy where esterification reactions occur at multiple stages of the synthesis.
    """
    # Count esterification reactions
    esterification_count = 0
    esterification_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal esterification_count, esterification_depths

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is an esterification reaction using predefined reaction types
                is_esterification = False

                # Check for known esterification reaction types
                if (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Transesterification", rsmi)
                    or checker.check_reaction("Oxidative esterification of primary alcohols", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                ):
                    is_esterification = True
                    print(f"Found predefined esterification reaction at depth {depth}")

                # If not found by reaction type, check by functional group transformation
                if not is_esterification:
                    try:
                        # Extract reactants and product
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant has a carboxylic acid group or an ester (for transesterification)
                        has_carboxylic_acid = False
                        has_reactant_ester = False
                        has_alcohol = False

                        for reactant in reactants:
                            if checker.check_fg("Carboxylic acid", reactant):
                                has_carboxylic_acid = True
                                print(f"Found carboxylic acid in reactant: {reactant}")
                            if checker.check_fg("Ester", reactant):
                                has_reactant_ester = True
                                print(f"Found ester in reactant: {reactant}")
                            if (
                                checker.check_fg("Primary alcohol", reactant)
                                or checker.check_fg("Secondary alcohol", reactant)
                                or checker.check_fg("Tertiary alcohol", reactant)
                            ):
                                has_alcohol = True
                                print(f"Found alcohol in reactant: {reactant}")

                        # Check if product has an ester group
                        has_ester = checker.check_fg("Ester", product)
                        if has_ester:
                            print(f"Found ester in product: {product}")

                        # If (reactant has carboxylic acid OR reactant has ester) and product has ester,
                        # it's likely an esterification or transesterification
                        if (
                            (has_carboxylic_acid and has_alcohol) or has_reactant_ester
                        ) and has_ester:
                            is_esterification = True
                            print(
                                f"Detected esterification by functional group change at depth {depth}"
                            )
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")

                if is_esterification:
                    esterification_count += 1
                    esterification_depths.append(depth)
                    print(f"Confirmed esterification reaction at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Total esterification reactions found: {esterification_count} at depths {esterification_depths}"
    )

    # Strategy is present if multiple esterification reactions are found
    return esterification_count >= 2
