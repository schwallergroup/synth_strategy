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
    Detects if the synthesis involves late-stage installation of a sulfonamide group.
    Looks for N-S bond formation between an amine and sulfonyl chloride in the late stages.
    """
    print(f"Analyzing route for late-stage sulfonamide formation")
    late_stage_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_sulfonamide

        if node["type"] == "reaction" and depth <= 3:  # Consider reactions within first four steps
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is a sulfonamide formation reaction - try specific reactions first
                is_primary_sulfonamide = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                )
                is_secondary_sulfonamide = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )

                print(
                    f"Reaction checks at depth {depth}: primary={is_primary_sulfonamide}, secondary={is_secondary_sulfonamide}"
                )

                # If specific reaction checks fail, look for general pattern
                is_sulfonamide_formation = is_primary_sulfonamide or is_secondary_sulfonamide

                # Verify sulfonyl chloride in reactants
                has_sulfonyl_chloride = any(
                    checker.check_fg("Sulfonyl halide", reactant) for reactant in reactants
                )
                if has_sulfonyl_chloride:
                    print(f"Found sulfonyl chloride in reactants at depth {depth}")

                # Verify amine in reactants
                has_amine = any(
                    checker.check_fg("Primary amine", reactant)
                    or checker.check_fg("Secondary amine", reactant)
                    for reactant in reactants
                )
                if has_amine:
                    print(f"Found amine in reactants at depth {depth}")

                # Verify sulfonamide in product
                has_sulfonamide_product = checker.check_fg("Sulfonamide", product)
                if has_sulfonamide_product:
                    print(f"Found sulfonamide in product at depth {depth}")

                # Verify sulfonamide was not already present in reactants
                sulfonamide_in_reactants = any(
                    checker.check_fg("Sulfonamide", reactant) for reactant in reactants
                )
                if sulfonamide_in_reactants:
                    print(f"Sulfonamide already present in reactants at depth {depth}")

                # If reaction check failed but we have all the right components, consider it a sulfonamide formation
                if (
                    not is_sulfonamide_formation
                    and has_sulfonyl_chloride
                    and has_amine
                    and has_sulfonamide_product
                    and not sulfonamide_in_reactants
                ):
                    print(f"Detected sulfonamide formation by component analysis at depth {depth}")
                    is_sulfonamide_formation = True

                # Confirm all conditions are met
                if is_sulfonamide_formation or (
                    has_sulfonyl_chloride
                    and has_amine
                    and has_sulfonamide_product
                    and not sulfonamide_in_reactants
                ):
                    print(f"Confirmed late-stage sulfonamide formation at depth {depth}")
                    late_stage_sulfonamide = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return late_stage_sulfonamide
