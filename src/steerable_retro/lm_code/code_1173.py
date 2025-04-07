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
    This function detects if the final step in the synthesis is a sulfonamide formation.
    """
    sulfonamide_at_depth_zero = False

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_at_depth_zero

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth == 1:  # Final reaction step is at depth 1
            # Extract reactants and product
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    print("No reaction SMILES found in metadata")
                    return

                print(f"Checking reaction at depth {depth}")
                print(f"RSMI: {rsmi}")

                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Reactants: {reactants}")
                print(f"Product: {product}")

                # Check if this is a sulfonamide formation reaction using reaction checker
                is_sulfonamide_reaction = (
                    checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    )
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                )

                print(f"Direct reaction check result: {is_sulfonamide_reaction}")

                # If direct reaction check fails, check for functional groups
                if not is_sulfonamide_reaction:
                    # Check for sulfonyl chloride in reactants
                    reactant_has_sulfonyl_chloride = any(
                        checker.check_fg("Sulfonyl halide", r) for r in reactants
                    )
                    print(f"Reactant has sulfonyl halide: {reactant_has_sulfonyl_chloride}")

                    # Check for primary or secondary amine in reactants
                    reactant_has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )
                    print(f"Reactant has amine: {reactant_has_amine}")

                    # Check for sulfonamide in product
                    product_has_sulfonamide = checker.check_fg("Sulfonamide", product)
                    print(f"Product has sulfonamide: {product_has_sulfonamide}")

                    # Determine if this is a sulfonamide formation
                    is_sulfonamide_reaction = (
                        reactant_has_sulfonyl_chloride
                        and reactant_has_amine
                        and product_has_sulfonamide
                    )

                if is_sulfonamide_reaction:
                    sulfonamide_at_depth_zero = True
                    print("Found sulfonamide formation at final step (depth 1)")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage sulfonamide formation: {sulfonamide_at_depth_zero}")
    return sulfonamide_at_depth_zero
