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
    This function detects if the synthesis route contains an oxidation
    of a primary alcohol to an aldehyde.
    """
    oxidation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal oxidation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Check if this is an oxidation reaction (forward direction)
            is_oxidation = (
                checker.check_reaction("Oxidation of alcohol to aldehyde", rsmi)
                or checker.check_reaction(
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                )
                or checker.check_reaction("Oxidation of alkene to aldehyde", rsmi)
            )

            if is_oxidation:
                print(f"Oxidation reaction detected: {rsmi}")

                # Check if reactants contain primary alcohol and product contains aldehyde
                has_primary_alcohol = any(
                    checker.check_fg("Primary alcohol", r) for r in reactants if r
                )
                has_aldehyde = checker.check_fg("Aldehyde", product) or checker.check_fg(
                    "Formaldehyde", product
                )

                if has_primary_alcohol and has_aldehyde:
                    print(f"Primary alcohol in reactants: {has_primary_alcohol}")
                    print(f"Aldehyde in product: {has_aldehyde}")
                    print("Alcohol oxidation to aldehyde confirmed")
                    oxidation_found = True

            # Check for the reverse reaction (reduction of aldehyde to alcohol)
            # Since we're traversing in retrosynthetic direction
            is_reduction = checker.check_reaction(
                "Reduction of aldehydes and ketones to alcohols", rsmi
            )

            if is_reduction:
                print(f"Reduction reaction detected: {rsmi}")

                has_aldehyde_reactant = any(
                    checker.check_fg("Aldehyde", r) or checker.check_fg("Formaldehyde", r)
                    for r in reactants
                    if r
                )
                has_primary_alcohol_product = checker.check_fg("Primary alcohol", product)

                if has_aldehyde_reactant and has_primary_alcohol_product:
                    print(f"Aldehyde in reactants: {has_aldehyde_reactant}")
                    print(f"Primary alcohol in product: {has_primary_alcohol_product}")
                    print("Aldehyde reduction to alcohol detected (reverse of target reaction)")
                    oxidation_found = True

            # Alternative approach: check for functional group changes even if reaction type not recognized
            if not oxidation_found and not is_oxidation and not is_reduction:
                # Forward direction: alcohol to aldehyde
                has_primary_alcohol_reactant = any(
                    checker.check_fg("Primary alcohol", r) for r in reactants if r
                )
                has_aldehyde_product = checker.check_fg("Aldehyde", product) or checker.check_fg(
                    "Formaldehyde", product
                )

                # Reverse direction: aldehyde to alcohol
                has_aldehyde_reactant = any(
                    checker.check_fg("Aldehyde", r) or checker.check_fg("Formaldehyde", r)
                    for r in reactants
                    if r
                )
                has_primary_alcohol_product = checker.check_fg("Primary alcohol", product)

                if (has_primary_alcohol_reactant and has_aldehyde_product) or (
                    has_aldehyde_reactant and has_primary_alcohol_product
                ):
                    print(
                        "Potential alcohol oxidation/reduction detected based on functional groups"
                    )

                    # Try to verify it's an oxidation/reduction by checking for absence of other reactions
                    # that might convert alcohol to aldehyde or vice versa
                    is_other_reaction = (
                        checker.check_reaction("Wohl-Ziegler bromination", rsmi)
                        or checker.check_reaction("Esterification", rsmi)
                        or checker.check_reaction(
                            "Reduction of carboxylic acid to primary alcohol", rsmi
                        )
                        or checker.check_reaction("Reduction of ester to primary alcohol", rsmi)
                        or checker.check_reaction("Aldehyde or ketone acetalization", rsmi)
                        or checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)
                    )

                    if not is_other_reaction:
                        print("Alcohol oxidation/reduction confirmed by functional group change")
                        oxidation_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return oxidation_found
