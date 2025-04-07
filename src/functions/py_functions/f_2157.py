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
    This function detects if a nitrile group is maintained throughout the synthesis.
    Returns True if:
    1. The final target molecule contains a nitrile group
    2. In each reaction where a nitrile is present in the product, it was also present in at least one reactant
    3. The nitrile group is maintained in at least 75% of the reactions
    """
    # Track reactions where nitrile is maintained
    reactions_with_nitrile_maintained = 0
    total_reactions_with_nitrile_product = 0

    # Check if the final target molecule contains a nitrile
    target_has_nitrile = checker.check_fg("Nitrile", route["smiles"])
    print(f"Target molecule has nitrile: {target_has_nitrile}")

    def dfs_traverse(node, depth=0):
        nonlocal reactions_with_nitrile_maintained, total_reactions_with_nitrile_product

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has nitrile
                product_has_nitrile = checker.check_fg("Nitrile", product_smiles)

                if product_has_nitrile:
                    total_reactions_with_nitrile_product += 1

                    # Check if any reactant has nitrile
                    reactants_have_nitrile = any(
                        checker.check_fg("Nitrile", reactant)
                        for reactant in reactants_smiles
                    )

                    # Nitrile is maintained if it's in both product and at least one reactant
                    if reactants_have_nitrile:
                        reactions_with_nitrile_maintained += 1
                        print(f"Nitrile group maintained in reaction: {rsmi}")
                    else:
                        print(f"Nitrile group created in reaction: {rsmi}")
                elif any(
                    checker.check_fg("Nitrile", reactant)
                    for reactant in reactants_smiles
                ):
                    print(f"Nitrile group lost in reaction: {rsmi}")

            except Exception as e:
                print(f"Error in nitrile detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True if:
    # 1. The target molecule has a nitrile group
    # 2. At least 75% of reactions maintain the nitrile group
    nitrile_maintenance_ratio = (
        reactions_with_nitrile_maintained / total_reactions_with_nitrile_product
        if total_reactions_with_nitrile_product > 0
        else 1.0
    )
    print(
        f"Nitrile maintenance ratio: {nitrile_maintenance_ratio} ({reactions_with_nitrile_maintained}/{total_reactions_with_nitrile_product})"
    )

    return target_has_nitrile and nitrile_maintenance_ratio >= 0.75
