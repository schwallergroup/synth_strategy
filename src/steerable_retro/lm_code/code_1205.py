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
    Detects if the route uses late-stage amine diversification as the final step.
    Specifically, looks for C-N bond formation in the last step of the synthesis.
    """
    found_late_amine = False
    final_reactions = []

    # First pass to identify final reactions (those leading directly to the target molecule)
    def identify_final_reactions(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Check if this reaction leads to a final product (the target molecule)
            for child in node.get("children", []):
                if child.get("type") == "mol" and not child.get("in_stock", False):
                    # If this child is the target molecule (has no further reactions)
                    has_reaction_child = False
                    for grandchild in child.get("children", []):
                        if grandchild.get("type") == "reaction":
                            has_reaction_child = True
                            break

                    if not has_reaction_child:
                        # This is a final reaction (leads to target molecule)
                        final_reactions.append((node, depth))
                        break

        # Continue traversal
        for child in node.get("children", []):
            identify_final_reactions(child, depth + 1)

    # Second pass to check final reactions for amine diversification
    def check_amine_diversification(reaction_node):
        nonlocal found_late_amine

        if "metadata" in reaction_node and "rsmi" in reaction_node["metadata"]:
            rsmi = reaction_node["metadata"]["rsmi"]
            print(f"Checking reaction: {rsmi}")

            # Check for amine diversification reactions
            amine_diversification_reactions = [
                "N-alkylation of primary amines with alkyl halides",
                "N-alkylation of secondary amines with alkyl halides",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                "Reductive amination with aldehyde",
                "Reductive amination with ketone",
                "Reductive amination with alcohol",
                "Alkylation of amines",
                "Buchwald-Hartwig",
                "{Buchwald-Hartwig}",
                "Negishi coupling",
                "Heck reaction with vinyl ester and amine",
                "Goldberg coupling",
                "Goldberg coupling aryl amine-aryl chloride",
                "Goldberg coupling aryl amide-aryl chloride",
                "Ullmann-Goldberg Substitution amine",
                "Ugi reaction",
                "Reductive amination",
                "{reductive amination}",
                "aza-Michael addition aromatic",
                "aza-Michael addition secondary",
                "aza-Michael addition primary",
                "Displacement of ethoxy group by primary amine",
                "Displacement of ethoxy group by secondary amine",
                "Eschweiler-Clarke Primary Amine Methylation",
                "Eschweiler-Clarke Secondary Amine Methylation",
                "Reductive methylation of primary amine with formaldehyde",
                "N-methylation",
                "Methylation with MeI_primary",
                "Methylation with MeI_secondary",
                "Methylation with MeI_tertiary",
                "DMS Amine methylation",
            ]

            for reaction_type in amine_diversification_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found amine diversification reaction: {reaction_type}")
                    found_late_amine = True
                    return

            # If no specific reaction type matched, check for reactants and products
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]
            reactants = reactants_part.split(".")

            # Check for amine functional groups in reactants
            has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
            has_secondary_amine = any(checker.check_fg("Secondary amine", r) for r in reactants)
            has_tertiary_amine_reactant = any(
                checker.check_fg("Tertiary amine", r) for r in reactants
            )

            # Check for leaving groups in reactants
            has_leaving_group = any(
                checker.check_fg("Primary halide", r)
                or checker.check_fg("Secondary halide", r)
                or checker.check_fg("Tertiary halide", r)
                or checker.check_fg("Aromatic halide", r)
                or checker.check_fg("Triflate", r)
                or checker.check_fg("Mesylate", r)
                or checker.check_fg("Tosylate", r)
                for r in reactants
            )

            # Check for amine groups in product
            has_primary_amine_product = checker.check_fg("Primary amine", product_part)
            has_secondary_amine_product = checker.check_fg("Secondary amine", product_part)
            has_tertiary_amine_product = checker.check_fg("Tertiary amine", product_part)

            # Check for C-N bond formation
            if (has_primary_amine or has_secondary_amine) and has_leaving_group:
                # Primary → Secondary or Secondary → Tertiary
                if (
                    has_primary_amine and has_secondary_amine_product and not has_secondary_amine
                ) or (
                    has_secondary_amine
                    and has_tertiary_amine_product
                    and not has_tertiary_amine_reactant
                ):
                    print("Found amine diversification: Amine alkylation pattern")
                    found_late_amine = True
                    return

            # Check for other amine formation patterns
            if not found_late_amine:
                # Check for amide coupling reactions
                if (
                    any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                    and any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )
                    and (
                        checker.check_fg("Primary amide", product_part)
                        or checker.check_fg("Secondary amide", product_part)
                        or checker.check_fg("Tertiary amide", product_part)
                    )
                ):
                    print("Found amine diversification: Amide coupling pattern")
                    found_late_amine = True
                    return

                # Check for urea formation
                if (
                    any(checker.check_fg("Isocyanate", r) for r in reactants)
                    and any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )
                    and checker.check_fg("Urea", product_part)
                ):
                    print("Found amine diversification: Urea formation pattern")
                    found_late_amine = True
                    return

                # Check for sulfonamide formation
                if (
                    any(checker.check_fg("Sulfonyl halide", r) for r in reactants)
                    and any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )
                    and checker.check_fg("Sulfonamide", product_part)
                ):
                    print("Found amine diversification: Sulfonamide formation pattern")
                    found_late_amine = True
                    return

    # Start traversal to identify final reactions
    identify_final_reactions(route)

    # Sort final reactions by depth (lowest depth = latest stage)
    final_reactions.sort(key=lambda x: x[1])

    # Check each final reaction for amine diversification, starting with the latest stage
    for reaction, _ in final_reactions:
        check_amine_diversification(reaction)
        if found_late_amine:
            break

    # If no final reactions were found or none had amine diversification,
    # try to find the last reaction in the route
    if not final_reactions or not found_late_amine:

        def find_last_reaction(node, current_depth=0, min_depth=float("inf"), last_reaction=None):
            if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
                if current_depth <= min_depth:
                    min_depth = current_depth
                    last_reaction = node

            for child in node.get("children", []):
                result_depth, result_reaction = find_last_reaction(
                    child, current_depth + 1, min_depth, last_reaction
                )
                if result_depth <= min_depth and result_reaction is not None:
                    min_depth = result_depth
                    last_reaction = result_reaction

            return min_depth, last_reaction

        _, last_reaction = find_last_reaction(route)
        if last_reaction:
            print("Checking last reaction in route")
            check_amine_diversification(last_reaction)

    return found_late_amine
