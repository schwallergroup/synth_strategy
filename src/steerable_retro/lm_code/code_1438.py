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
    This function detects if the synthesis involves oxidation of a primary alcohol
    to an aldehyde.
    """
    found_oxidation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_oxidation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            # Check if this is an oxidation reaction
            if checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi):
                print(f"Found oxidation reaction at depth {depth}: {rsmi}")
                found_oxidation = True
                return

            # Alternative check using "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones"
            if checker.check_reaction(
                "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
            ):
                print(f"Found alcohol oxidation reaction at depth {depth}: {rsmi}")

                # Check if reactants contain primary alcohol and product contains aldehyde
                reactants = reactants_str.split(".")

                has_primary_alcohol = False
                for reactant in reactants:
                    if checker.check_fg("Primary alcohol", reactant):
                        has_primary_alcohol = True
                        break

                has_aldehyde = checker.check_fg("Aldehyde", product_str)

                if has_primary_alcohol and has_aldehyde:
                    print(f"Confirmed primary alcohol to aldehyde oxidation at depth {depth}")
                    found_oxidation = True
                    return

            # Fallback method: check for functional group transformation directly
            try:
                reactants = reactants_str.split(".")

                # Check for primary alcohol in reactants
                primary_alcohol_found = False
                for reactant in reactants:
                    if checker.check_fg("Primary alcohol", reactant):
                        primary_alcohol_found = True
                        break

                # Check for aldehyde in product
                aldehyde_found = checker.check_fg("Aldehyde", product_str)

                if primary_alcohol_found and aldehyde_found:
                    # Additional check to ensure it's an oxidation reaction
                    # This is a simplified check - in a real scenario, we would verify atom mapping
                    print(f"Detected primary alcohol to aldehyde transformation at depth {depth}")
                    found_oxidation = True
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_oxidation
