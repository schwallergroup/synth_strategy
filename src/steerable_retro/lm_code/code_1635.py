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
    Detects a synthetic strategy involving late-stage formation of a sulfonamide group.
    Late-stage means in the second half of the synthesis (higher depth values in retrosynthesis).
    """
    found_sulfonamide_formation = False
    sulfonamide_formation_depth = float("inf")  # Initialize to infinity
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_sulfonamide_formation, sulfonamide_formation_depth, max_depth

        # Track maximum depth to determine synthesis length
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")

            # Check if product has sulfonamide
            product_has_sulfonamide = checker.check_fg("Sulfonamide", product_part)

            # Check if this is a sulfonamide formation reaction
            is_sulfonamide_reaction = (
                checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                )
                or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )
                or checker.check_reaction("sulfon_amide", rsmi)
                or checker.check_reaction("Schotten-Baumann to ester", rsmi)
            )

            if is_sulfonamide_reaction and product_has_sulfonamide:
                # Check if any reactant doesn't have the sulfonamide
                has_reactant_without_sulfonamide = any(
                    not checker.check_fg("Sulfonamide", reactant) for reactant in reactants
                )

                if has_reactant_without_sulfonamide:
                    # If we find a sulfonamide formation at a lower depth (later stage),
                    # or if this is the first one we've found, record it
                    if not found_sulfonamide_formation or depth < sulfonamide_formation_depth:
                        found_sulfonamide_formation = True
                        sulfonamide_formation_depth = depth
                        print(f"Found sulfonamide formation at depth {depth}, reaction: {rsmi}")

            # Additional check for sulfonamide formation by looking at functional groups
            elif product_has_sulfonamide:
                # Check if any reactant doesn't have the sulfonamide
                has_reactant_without_sulfonamide = any(
                    not checker.check_fg("Sulfonamide", reactant) for reactant in reactants
                )

                if has_reactant_without_sulfonamide:
                    # Check if reactants have sulfonyl chloride and amine
                    has_sulfonyl_chloride = any(
                        checker.check_fg("Sulfonyl halide", reactant) for reactant in reactants
                    )
                    has_amine = any(
                        checker.check_fg("Primary amine", reactant)
                        or checker.check_fg("Secondary amine", reactant)
                        or checker.check_fg("Aniline", reactant)
                        for reactant in reactants
                    )

                    if has_sulfonyl_chloride and has_amine:
                        # This is likely a sulfonamide formation reaction not captured by the reaction checks
                        if not found_sulfonamide_formation or depth < sulfonamide_formation_depth:
                            found_sulfonamide_formation = True
                            sulfonamide_formation_depth = depth
                            print(
                                f"Found sulfonamide formation via functional groups at depth {depth}"
                            )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if sulfonamide formation is late-stage (in the second half of synthesis)
    # In retrosynthesis, higher depth values are later in the forward synthesis
    late_stage = False
    if found_sulfonamide_formation and max_depth > 0:
        midpoint = max_depth / 2
        # Late stage means depth is greater than or equal to midpoint (deeper in the tree)
        late_stage = sulfonamide_formation_depth >= midpoint
        print(
            f"Max depth: {max_depth}, Midpoint: {midpoint}, Sulfonamide formation depth: {sulfonamide_formation_depth}"
        )
        print(
            f"Is late stage? {late_stage} (depth {sulfonamide_formation_depth} >= midpoint {midpoint})"
        )
    else:
        print(f"Max depth: {max_depth}, Sulfonamide formation not found or max_depth is 0")

    return late_stage
