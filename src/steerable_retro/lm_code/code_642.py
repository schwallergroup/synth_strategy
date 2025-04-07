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
    This function detects if nitrogen functionalization (reductive amination, amide formation,
    N-alkylation, etc.) occurs in the late stage of the synthesis (last 40% of steps).
    """
    total_steps = 0
    n_functionalization_depths = []

    # First count total steps
    def count_steps(node):
        nonlocal total_steps
        if node["type"] == "reaction":
            total_steps += 1
        for child in node.get("children", []):
            count_steps(child)

    count_steps(route)

    if total_steps == 0:
        print("No reaction steps found in the route")
        return False

    # Then identify nitrogen functionalization steps and their depths
    def find_n_functionalization(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check for various nitrogen functionalization reactions using the checker
            n_functionalization_reactions = [
                "Reductive amination with aldehyde",
                "Reductive amination with ketone",
                "Reductive amination with alcohol",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "N-alkylation of primary amines with alkyl halides",
                "N-alkylation of secondary amines with alkyl halides",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Carboxylic acid with primary amine to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Schotten-Baumann_amide",
                "Urea synthesis via isocyanate and primary amine",
                "Urea synthesis via isocyanate and secondary amine",
                "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
            ]

            for reaction_type in n_functionalization_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    n_functionalization_depths.append(depth)
                    print(f"Nitrogen functionalization found at depth {depth}: {reaction_type}")
                    print(f"Reaction SMILES: {rsmi}")
                    break

            # If no specific reaction type matched, check for general patterns
            if not any(checker.check_reaction(rxn, rsmi) for rxn in n_functionalization_reactions):
                # Extract reactants and product
                try:
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check if any reactant has nitrogen-containing functional groups
                    reactant_n_groups = [
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Amide",
                        "Nitrile",
                        "Nitro group",
                        "Azide",
                    ]

                    # Check if product has different or modified nitrogen functional groups
                    product_n_groups = [
                        "Primary amine",
                        "Secondary amine",
                        "Tertiary amine",
                        "Primary amide",
                        "Secondary amide",
                        "Tertiary amide",
                    ]

                    reactant_has_n = any(
                        any(checker.check_fg(fg, r) for fg in reactant_n_groups) for r in reactants
                    )
                    product_has_n = any(checker.check_fg(fg, product) for fg in product_n_groups)

                    # If both reactants and products have nitrogen groups, it might be a nitrogen functionalization
                    if reactant_has_n and product_has_n:
                        # Check if the nitrogen environment changed
                        reactant_n_types = set()
                        for r in reactants:
                            for fg in reactant_n_groups:
                                if checker.check_fg(fg, r):
                                    reactant_n_types.add(fg)

                        product_n_types = set()
                        for fg in product_n_groups:
                            if checker.check_fg(fg, product):
                                product_n_types.add(fg)

                        # If the nitrogen functional group types changed, it's likely a nitrogen functionalization
                        if reactant_n_types != product_n_types:
                            n_functionalization_depths.append(depth)
                            print(f"Potential nitrogen functionalization found at depth {depth}")
                            print(f"Reaction SMILES: {rsmi}")
                            print(f"Reactant N groups: {reactant_n_types}")
                            print(f"Product N groups: {product_n_types}")
                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

        for child in node.get("children", []):
            find_n_functionalization(child, depth + 1)

    find_n_functionalization(route)

    # Check if any nitrogen functionalization occurs in the last 40% of steps
    if not n_functionalization_depths:
        print("No nitrogen functionalization reactions found")
        return False

    late_stage_threshold = total_steps * 0.4
    for depth in n_functionalization_depths:
        if depth < late_stage_threshold:  # Lower depth values are later in the synthesis
            print(
                f"Late-stage nitrogen functionalization detected (depth {depth} out of {total_steps} steps)"
            )
            return True

    print(
        f"Nitrogen functionalizations found at depths {n_functionalization_depths}, but none in late stage (threshold: {late_stage_threshold})"
    )
    return False
