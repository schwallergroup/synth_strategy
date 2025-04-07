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
    Detects if the synthesis route involves oxidation of a primary alcohol to an aldehyde.
    """
    oxidation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal oxidation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reactants and product from the forward reaction
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            print(f"Depth {depth} - Examining reaction: {rsmi}")

            # Check if this is an oxidation reaction (forward direction)
            is_oxidation_reaction = checker.check_reaction(
                "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
            )

            if is_oxidation_reaction:
                print(f"  Found oxidation reaction: {rsmi}")
                # In forward direction: alcohol (reactant) -> aldehyde (product)
                has_primary_alcohol = any(
                    checker.check_fg("Primary alcohol", r) for r in reactants_smiles
                )
                has_aldehyde = checker.check_fg("Aldehyde", product_smiles)

                print(f"  Has primary alcohol in reactants: {has_primary_alcohol}")
                print(f"  Has aldehyde in product: {has_aldehyde}")

                if has_primary_alcohol and has_aldehyde:
                    print(f"  ✓ Confirmed alcohol oxidation to aldehyde: {rsmi}")
                    oxidation_found = True

            # Also check for reduction reaction (forward direction)
            # This would be an oxidation when considered in reverse
            is_reduction_reaction = checker.check_reaction(
                "Reduction of aldehydes and ketones to alcohols", rsmi
            )

            # Manual check for aldehyde to alcohol transformation
            if not is_reduction_reaction:
                for reactant in reactants_smiles:
                    if checker.check_fg("Aldehyde", reactant) and checker.check_fg(
                        "Primary alcohol", product_smiles
                    ):
                        print(f"  Detected aldehyde to alcohol transformation manually: {rsmi}")
                        is_reduction_reaction = True
                        break

            if is_reduction_reaction:
                print(f"  Found reduction reaction: {rsmi}")
                # In forward direction: aldehyde (reactant) -> alcohol (product)
                has_aldehyde_in_reactants = any(
                    checker.check_fg("Aldehyde", r) for r in reactants_smiles
                )
                has_primary_alcohol_in_product = checker.check_fg("Primary alcohol", product_smiles)

                print(f"  Has aldehyde in reactants: {has_aldehyde_in_reactants}")
                print(f"  Has primary alcohol in product: {has_primary_alcohol_in_product}")

                # Print which reactant has the aldehyde for debugging
                if has_aldehyde_in_reactants:
                    for i, r in enumerate(reactants_smiles):
                        if checker.check_fg("Aldehyde", r):
                            print(f"    Aldehyde found in reactant {i}: {r}")

                if has_aldehyde_in_reactants and has_primary_alcohol_in_product:
                    print(
                        f"  ✓ Confirmed aldehyde reduction to alcohol (reverse of oxidation): {rsmi}"
                    )
                    oxidation_found = True

            # Check for other potential oxidation reactions
            if not is_oxidation_reaction and not is_reduction_reaction:
                # Check for alcohol oxidation to carboxylic acid (might go through aldehyde)
                is_alcohol_to_acid = checker.check_reaction(
                    "Oxidation of alcohol to carboxylic acid", rsmi
                )

                if is_alcohol_to_acid:
                    print(f"  Found alcohol to acid oxidation (might involve aldehyde): {rsmi}")
                    has_primary_alcohol = any(
                        checker.check_fg("Primary alcohol", r) for r in reactants_smiles
                    )
                    has_carboxylic_acid = checker.check_fg("Carboxylic acid", product_smiles)

                    if has_primary_alcohol and has_carboxylic_acid:
                        # This likely went through an aldehyde intermediate
                        print(
                            f"  ✓ Confirmed alcohol oxidation to carboxylic acid (via aldehyde): {rsmi}"
                        )
                        oxidation_found = True

            # Direct pattern check for C=O to CH2-OH transformation or vice versa
            # This is a fallback in case the reaction checks miss something
            if not oxidation_found:
                # Look for pattern-based evidence of aldehyde-alcohol interconversion
                for reactant in reactants_smiles:
                    # Check for aldehyde pattern in atom-mapped SMILES
                    if "[CH]=[O]" in reactant and "[CH2][OH]" in product_smiles:
                        print(f"  ✓ Pattern-based detection of aldehyde reduction: {rsmi}")
                        oxidation_found = True
                        break
                    elif "[CH2][OH]" in reactant and "[CH]=[O]" in product_smiles:
                        print(f"  ✓ Pattern-based detection of alcohol oxidation: {rsmi}")
                        oxidation_found = True
                        break

        # Traverse children (retrosynthetic direction)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Alcohol oxidation to aldehyde strategy: {oxidation_found}")
    return oxidation_found
