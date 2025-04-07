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
    Detects a synthesis strategy with:
    1. Late-stage amide coupling (final step)
    2. Nitro reduction in the synthesis path
    3. Linear synthesis with convergent coupling only at the final step
    """
    has_amide_coupling_final_step = False
    has_nitro_reduction = False

    # Track the final step depth
    final_step_depth = None

    # First pass to identify the final step depth
    def find_final_step(node, depth=0):
        nonlocal final_step_depth

        if node["type"] == "reaction":
            if final_step_depth is None or depth < final_step_depth:
                final_step_depth = depth

        for child in node.get("children", []):
            find_final_step(child, depth + 1)

    # Find the final step depth
    find_final_step(route)
    print(f"Final step identified at depth: {final_step_depth}")

    # Second pass to analyze the synthesis
    def dfs_traverse(node, depth=0, path=None):
        nonlocal has_amide_coupling_final_step, has_nitro_reduction

        if path is None:
            path = []

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for amide coupling at the final step
                if depth == final_step_depth:
                    # Check for amide coupling reactions
                    amide_reaction_types = [
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Carboxylic acid with primary amine to amide",
                        "Ester with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        "Schotten-Baumann to ester",
                        "{Schotten-Baumann_amide}",
                        "Acyl chloride with secondary amine to amide",
                        "Ester with secondary amine to amide",
                        "Acyl chloride with ammonia to amide",
                    ]

                    is_amide_reaction = any(
                        checker.check_reaction(rxn, rsmi) for rxn in amide_reaction_types
                    )

                    # Verify reactants have carboxylic acid/derivative and amine
                    has_acid_derivative = False
                    has_amine = False

                    for reactant in reactants:
                        if (
                            checker.check_fg("Carboxylic acid", reactant)
                            or checker.check_fg("Acyl halide", reactant)
                            or checker.check_fg("Ester", reactant)
                            or checker.check_fg("Anhydride", reactant)
                        ):
                            has_acid_derivative = True
                            print(f"Found acid derivative in reactant: {reactant}")

                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        ):
                            has_amine = True
                            print(f"Found amine in reactant: {reactant}")

                    # Verify product has amide
                    has_amide_product = (
                        checker.check_fg("Primary amide", product_part)
                        or checker.check_fg("Secondary amide", product_part)
                        or checker.check_fg("Tertiary amide", product_part)
                    )

                    # If no specific reaction type is detected, check for functional group changes
                    if (
                        not is_amide_reaction
                        and has_acid_derivative
                        and has_amine
                        and has_amide_product
                    ):
                        is_amide_reaction = True
                        print("Detected amide formation based on functional group changes")

                    if is_amide_reaction and has_amide_product:
                        has_amide_coupling_final_step = True
                        print("Found amide coupling in final step")

                # Check for nitro reduction at any depth
                has_nitro_reactant = False
                for reactant in reactants:
                    if checker.check_fg("Nitro group", reactant):
                        has_nitro_reactant = True
                        print(f"Found nitro group in reactant: {reactant}")

                if has_nitro_reactant:
                    # Check if product has amine but no nitro group
                    has_amine_product = checker.check_fg(
                        "Primary amine", product_part
                    ) or checker.check_fg("Aniline", product_part)
                    has_nitro_product = checker.check_fg("Nitro group", product_part)

                    # Check for nitro reduction reaction or infer from functional group changes
                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        has_nitro_reduction = True
                        print(f"Found nitro reduction reaction at depth {depth}")
                    elif has_amine_product and not has_nitro_product and has_nitro_reactant:
                        has_nitro_reduction = True
                        print(
                            f"Inferred nitro reduction at depth {depth} based on functional group changes"
                        )

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Track this node in the current path
        current_path = path + [node]

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1, current_path)

    # Start traversal from root
    dfs_traverse(route)

    print(
        f"Final results - Amide coupling: {has_amide_coupling_final_step}, Nitro reduction: {has_nitro_reduction}"
    )

    # Return True if both conditions are met
    return has_amide_coupling_final_step and has_nitro_reduction
