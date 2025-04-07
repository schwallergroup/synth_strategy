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
    This function detects a synthetic strategy involving late-stage sulfonamide formation
    (typically in the final step of the synthesis).
    """
    found_sulfonamide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_sulfonamide_formation

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check if it's a late-stage reaction (depth 0 or 1)
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                reactants = reactants_part.split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a sulfonamide synthesis reaction
                if checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                ):
                    print(f"Found sulfonamide synthesis reaction: {rsmi}")
                    found_sulfonamide_formation = True
                else:
                    # Fallback to checking functional groups if reaction check fails
                    has_sulfonyl_chloride = False
                    has_amine = False
                    has_sulfonamide_product = False

                    # Check reactants for required functional groups
                    for reactant in reactants:
                        if checker.check_fg("Sulfonyl halide", reactant):
                            has_sulfonyl_chloride = True
                            print(f"Found sulfonyl halide in reactant: {reactant}")
                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            has_amine = True
                            print(f"Found amine in reactant: {reactant}")

                    # Check product for sulfonamide group
                    if checker.check_fg("Sulfonamide", product):
                        has_sulfonamide_product = True
                        print(f"Found sulfonamide in product: {product}")

                    # If all conditions are met, mark as found
                    if has_sulfonyl_chloride and has_amine and has_sulfonamide_product:
                        found_sulfonamide_formation = True
                        print(
                            f"Found late-stage sulfonamide formation through functional group analysis: {rsmi}"
                        )
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage sulfonamide formation strategy detected: {found_sulfonamide_formation}")
    return found_sulfonamide_formation
