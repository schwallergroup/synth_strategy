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
    Detects if the synthesis route uses a late-stage amide coupling strategy.
    This checks if an amide bond is formed in the final step (depth 0).
    """
    amide_formation_at_late_stage = False

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_at_late_stage

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Check depth 0 and 1 for late-stage
            print(f"Examining reaction at depth {depth}")

            # Get reaction SMILES
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                print("No reaction SMILES found")
                return

            print(f"Reaction SMILES: {rsmi}")

            # Extract reactants and product
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is an amide formation reaction using the checker function
            is_amide_reaction = False

            # Check for various amide formation reaction types
            amide_reaction_types = [
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Acyl chloride with ammonia to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Carboxylic acid with primary amine to amide",
                "Ester with ammonia to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Schotten-Baumann_amide",
            ]

            for reaction_type in amide_reaction_types:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Detected amide formation reaction: {reaction_type}")
                    is_amide_reaction = True
                    break

            # If no specific reaction type matched, check for functional group changes
            if not is_amide_reaction:
                # Check for carboxylic acid, acyl halide, or ester in reactants
                has_acid = any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants
                )
                has_acyl_halide = any(
                    checker.check_fg("Acyl halide", r) for r in reactants
                )
                has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants)

                # Check for amine in reactants
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants
                )
                has_secondary_amine = any(
                    checker.check_fg("Secondary amine", r) for r in reactants
                )

                # Check for amide in product
                has_amide_product = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                print(
                    f"Carboxylic acid: {has_acid}, Acyl halide: {has_acyl_halide}, Ester: {has_ester}, Anhydride: {has_anhydride}"
                )
                print(
                    f"Primary amine: {has_primary_amine}, Secondary amine: {has_secondary_amine}"
                )
                print(f"Amide in product: {has_amide_product}")

                # Check if we have the right combination of reactants and product
                if has_amide_product and (
                    (has_acid and (has_primary_amine or has_secondary_amine))
                    or (has_acyl_halide and (has_primary_amine or has_secondary_amine))
                    or (has_ester and (has_primary_amine or has_secondary_amine))
                    or (has_anhydride and (has_primary_amine or has_secondary_amine))
                ):
                    print(
                        f"Detected amide formation at depth {depth} through functional group analysis"
                    )
                    is_amide_reaction = True

            # If we found an amide formation reaction at the appropriate depth
            if is_amide_reaction:
                if depth == 0:
                    print("Confirmed late-stage amide coupling at final step (depth 0)")
                    amide_formation_at_late_stage = True
                elif (
                    depth == 1 and not amide_formation_at_late_stage
                ):  # Only set if depth 0 hasn't already set it
                    print(
                        "Confirmed late-stage amide coupling at penultimate step (depth 1)"
                    )
                    amide_formation_at_late_stage = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {amide_formation_at_late_stage}")
    return amide_formation_at_late_stage
