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
    This function detects if the synthetic route contains multiple alcohol oxidation steps.
    """
    oxidation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal oxidation_count

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for alcohol oxidation reactions directly
                if checker.check_reaction(
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                ) or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi):
                    oxidation_count += 1
                    print(f"Found alcohol oxidation reaction at depth {depth}: {rsmi}")
                else:
                    # In retrosynthesis, the product is the starting material and reactants are the target compounds
                    # So we need to check if the product has alcohol and reactants have oxidized forms

                    # Check if product has alcohol group
                    has_alcohol_in_product = (
                        checker.check_fg("Primary alcohol", product_smiles)
                        or checker.check_fg("Secondary alcohol", product_smiles)
                        or checker.check_fg("Tertiary alcohol", product_smiles)
                        or checker.check_fg("Aromatic alcohol", product_smiles)
                    )

                    if has_alcohol_in_product:
                        # Check if any reactant has aldehyde, ketone, or carboxylic acid
                        for reactant in reactants_smiles:
                            has_oxidized_group = (
                                checker.check_fg("Aldehyde", reactant)
                                or checker.check_fg("Formaldehyde", reactant)
                                or checker.check_fg("Ketone", reactant)
                                or checker.check_fg("Carboxylic acid", reactant)
                            )

                            if has_oxidized_group:
                                oxidation_count += 1
                                print(
                                    f"Found alcohol oxidation (by FG check) at depth {depth}: {rsmi}"
                                )
                                break

                    # Also check the forward direction for completeness
                    for reactant in reactants_smiles:
                        # Check if reactant has alcohol group
                        has_alcohol = (
                            checker.check_fg("Primary alcohol", reactant)
                            or checker.check_fg("Secondary alcohol", reactant)
                            or checker.check_fg("Tertiary alcohol", reactant)
                            or checker.check_fg("Aromatic alcohol", reactant)
                        )

                        if has_alcohol:
                            # Check if product has aldehyde, ketone, or carboxylic acid
                            has_oxidized_product = (
                                checker.check_fg("Aldehyde", product_smiles)
                                or checker.check_fg("Formaldehyde", product_smiles)
                                or checker.check_fg("Ketone", product_smiles)
                                or checker.check_fg("Carboxylic acid", product_smiles)
                            )

                            if has_oxidized_product:
                                oxidation_count += 1
                                print(
                                    f"Found alcohol oxidation (forward direction) at depth {depth}: {rsmi}"
                                )
                                break
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple oxidation steps are found
    result = oxidation_count >= 2
    print(f"Multiple alcohol oxidation strategy detected: {result} (count: {oxidation_count})")
    return result
