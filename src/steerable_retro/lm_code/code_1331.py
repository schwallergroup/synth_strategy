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
    This function detects an alcohol oxidation strategy where:
    A secondary alcohol is oxidized to a ketone
    """
    alcohol_oxidation_found = False

    def dfs_traverse(node, depth=0):
        nonlocal alcohol_oxidation_found

        print(f"Traversing node at depth {depth}: {node.get('type', 'unknown')}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            # In retrosynthesis, the product is on the left and reactants on the right
            product = rsmi.split(">")[0]
            reactants = rsmi.split(">")[-1].split(".")

            print(f"Checking reaction: {rsmi}")

            # Check if this is an oxidation reaction
            is_oxidation = checker.check_reaction(
                "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
            )

            # Also check for other oxidation reactions that might convert secondary alcohols to ketones
            if not is_oxidation:
                # Check if the reaction involves a secondary alcohol being oxidized to a ketone
                has_secondary_alcohol = False
                has_ketone = False

                # Check if product has a ketone
                has_ketone = checker.check_fg("Ketone", product) if product else False
                if has_ketone:
                    print(f"Found product with ketone: {product}")

                # Check if any reactant has a secondary alcohol
                for r in reactants:
                    if not r:
                        continue

                    if checker.check_fg("Secondary alcohol", r):
                        has_secondary_alcohol = True
                        print(f"Found reactant with secondary alcohol: {r}")
                        break

                if has_secondary_alcohol and has_ketone:
                    # This is likely a secondary alcohol oxidation to ketone
                    is_oxidation = True

            if is_oxidation:
                print(f"Found potential alcohol oxidation reaction: {rsmi}")

                # Check for secondary alcohol in reactants
                has_secondary_alcohol = False

                for r in reactants:
                    if not r:
                        continue

                    if checker.check_fg("Secondary alcohol", r):
                        has_secondary_alcohol = True
                        print(f"Found reactant with secondary alcohol: {r}")

                # Check for ketone in product
                has_ketone = checker.check_fg("Ketone", product) if product else False

                if has_ketone:
                    print(f"Found product with ketone: {product}")

                # Confirm secondary alcohol to ketone transformation
                if has_secondary_alcohol and has_ketone:
                    alcohol_oxidation_found = True
                    print("Confirmed alcohol oxidation strategy: secondary alcohol to ketone")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return alcohol_oxidation_found
