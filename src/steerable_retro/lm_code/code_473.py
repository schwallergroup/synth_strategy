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
    Detects if the synthesis route includes a late-stage N-alkylation
    (in the second half of the synthesis depth).
    """
    max_depth = 0
    n_alkylation_depths = []

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to find N-alkylation
    def find_n_alkylation(node, current_depth=0):
        nonlocal n_alkylation_depths

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for N-alkylation reaction types
            n_alkylation_reactions = [
                "N-alkylation of primary amines with alkyl halides",
                "N-alkylation of secondary amines with alkyl halides",
                "Alkylation of amines",
                "Methylation with MeI_primary",
                "Methylation with MeI_secondary",
                "Methylation with MeI_tertiary",
                "DMS Amine methylation",
                "N-methylation",
                "Eschweiler-Clarke Primary Amine Methylation",
                "Eschweiler-Clarke Secondary Amine Methylation",
                "Reductive methylation of primary amine with formaldehyde",
                "Reductive amination with aldehyde",
                "Reductive amination with ketone",
                "Reductive amination with alcohol",
                "Mitsunobu_sulfonamide",
                "Mitsunobu_tetrazole_1",
                "Mitsunobu_tetrazole_2",
                "Mitsunobu_tetrazole_3",
                "Mitsunobu_tetrazole_4",
            ]

            is_n_alkylation = False
            for rxn_type in n_alkylation_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(
                        f"Found N-alkylation reaction '{rxn_type}' at depth {current_depth} (late-stage depth: {max_depth - current_depth})"
                    )
                    is_n_alkylation = True
                    break

            # If no specific reaction type matched, try to detect N-alkylation pattern
            if not is_n_alkylation:
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for primary/secondary amine in reactants
                    reactants_have_amine = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary amine", reactant)
                            or checker.check_fg("Secondary amine", reactant)
                            or checker.check_fg("Aniline", reactant)
                        ):
                            reactants_have_amine = True
                            break

                    # Check for more substituted amine in product
                    product_has_alkylated_amine = False
                    if (
                        checker.check_fg("Secondary amine", product)
                        or checker.check_fg("Tertiary amine", product)
                        or checker.check_fg("Sulfonamide", product)
                    ):
                        product_has_alkylated_amine = True

                    # Additional check: ensure there's an alkyl halide or similar alkylating agent
                    reactants_have_alkylating_agent = False
                    for reactant in reactants:
                        if (
                            checker.check_fg("Primary halide", reactant)
                            or checker.check_fg("Secondary halide", reactant)
                            or checker.check_fg("Tertiary halide", reactant)
                            or checker.check_fg("Tosylate", reactant)
                            or checker.check_fg("Mesylate", reactant)
                            or checker.check_fg("Triflate", reactant)
                            or checker.check_fg("Aldehyde", reactant)  # For reductive amination
                            or checker.check_fg("Ketone", reactant)
                        ):  # For reductive amination
                            reactants_have_alkylating_agent = True
                            break

                    if (
                        reactants_have_amine
                        and product_has_alkylated_amine
                        and reactants_have_alkylating_agent
                    ):
                        print(
                            f"Found N-alkylation pattern at depth {current_depth} (late-stage depth: {max_depth - current_depth})"
                        )
                        is_n_alkylation = True
                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

            if is_n_alkylation:
                # In retrosynthesis, lower depths are later stages
                late_stage_depth = max_depth - current_depth
                n_alkylation_depths.append(late_stage_depth)

        for child in node.get("children", []):
            find_n_alkylation(child, current_depth + 1)

    find_max_depth(route)
    print(f"Maximum synthesis depth: {max_depth}")
    find_n_alkylation(route)

    # Check if any N-alkylation occurs in the second half of the synthesis (late-stage)
    if n_alkylation_depths:
        # Find the latest N-alkylation (highest late-stage depth)
        latest_n_alkylation = max(n_alkylation_depths)
        ratio = latest_n_alkylation / max_depth if max_depth > 0 else 0
        print(
            f"Latest N-alkylation at late-stage depth {latest_n_alkylation} out of {max_depth}, ratio: {ratio:.2f}"
        )

        # Late stage is defined as occurring in the second half of the synthesis
        if ratio >= 0.5:
            print(f"N-alkylation is late-stage (depth ratio: {ratio:.2f})")
            return True
        else:
            print(f"N-alkylation is early-stage (depth ratio: {ratio:.2f})")
            return False
    else:
        print("No N-alkylation reactions found in the synthesis route")
        return False
