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
    This function detects late-stage sulfonamide formation in the synthetic route.
    """
    sulfonamide_formation_detected = False
    depth_of_formation = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formation_detected, depth_of_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Try to detect sulfonamide formation
            try:
                # Check if product contains sulfonamide group
                if checker.check_fg("Sulfonamide", product):
                    print(f"Product contains sulfonamide at depth {depth}")

                    # Check if any reactant already contains sulfonamide
                    reactants_with_sulfonamide = [
                        r for r in reactants if checker.check_fg("Sulfonamide", r)
                    ]

                    if not reactants_with_sulfonamide:
                        print("No reactants contain sulfonamide - potential formation reaction")

                        # Check if this is a known sulfonamide formation reaction
                        is_sulfonamide_reaction = (
                            checker.check_reaction(
                                "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                            )
                            or checker.check_reaction(
                                "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                            )
                            or checker.check_reaction(
                                "Schotten-Baumann to ester", rsmi
                            )  # Sometimes mislabeled
                        )

                        # Verify reactants: need sulfonyl chloride and amine
                        has_sulfonyl_chloride = any(
                            checker.check_fg("Sulfonyl halide", r) for r in reactants
                        )
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants
                        )

                        print(
                            f"Reaction check: {is_sulfonamide_reaction}, Sulfonyl halide: {has_sulfonyl_chloride}, Amine: {has_amine}"
                        )

                        if is_sulfonamide_reaction or (has_sulfonyl_chloride and has_amine):
                            print(f"Sulfonamide formation confirmed at depth {depth}")
                            sulfonamide_formation_detected = True
                            # Track the depth to determine if it's late-stage (lower depth = later stage)
                            depth_of_formation = min(depth_of_formation, depth)
            except Exception as e:
                print(f"Error in sulfonamide detection: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(
        f"Final result: formation detected: {sulfonamide_formation_detected}, at depth: {depth_of_formation}"
    )
    # Consider it late-stage if it occurs in the first few steps of the synthesis
    # Increased threshold to depth 4 to be more inclusive
    return sulfonamide_formation_detected and depth_of_formation <= 4
