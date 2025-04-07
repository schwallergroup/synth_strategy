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
    This function detects a synthesis strategy that involves multiple reduction
    and oxidation steps in sequence.
    """
    # Track reaction steps with their types and depths
    reaction_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]

                # Check for reduction reactions using reaction checkers
                is_reduction = False
                if (
                    checker.check_reaction("Reduction of ester to primary alcohol", rsmi)
                    or checker.check_reaction(
                        "Reduction of carboxylic acid to primary alcohol", rsmi
                    )
                    or checker.check_reaction("Reduction of nitrile to amine", rsmi)
                    or checker.check_reaction(
                        "Reduction of aldehydes and ketones to alcohols", rsmi
                    )
                    or checker.check_reaction("Reduction of primary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of secondary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of tertiary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                ):
                    is_reduction = True
                    reaction_sequence.append(("reduction", depth))
                    print(f"Detected reduction reaction at depth {depth}: {rsmi}")

                # Check for oxidation reactions using reaction checkers
                is_oxidation = False
                if (
                    checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                    or checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                    )
                    or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of ketone to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of alkene to aldehyde", rsmi)
                    or checker.check_reaction("Oxidation of boronic acids", rsmi)
                    or checker.check_reaction("Oxidation of boronic esters", rsmi)
                ):
                    is_oxidation = True
                    reaction_sequence.append(("oxidation", depth))
                    print(f"Detected oxidation reaction at depth {depth}: {rsmi}")

                # If no specific reaction type was detected, check for functional group changes
                if not is_reduction and not is_oxidation:
                    reactants_part = rsmi.split(">")[0]
                    product_part = rsmi.split(">")[-1]

                    # Get reactants and product SMILES
                    reactants_smiles = reactants_part.split(".")
                    product_smiles = product_part

                    # Additional reduction checks based on functional groups
                    if any(
                        checker.check_fg("Nitrile", r) for r in reactants_smiles
                    ) and checker.check_fg("Primary amine", product_smiles):
                        reaction_sequence.append(("reduction", depth))
                        print(f"Detected reduction: Nitrile to amine at depth {depth}")

                    elif any(checker.check_fg("Ester", r) for r in reactants_smiles) and (
                        checker.check_fg("Primary alcohol", product_smiles)
                        or checker.check_fg("Secondary alcohol", product_smiles)
                        or checker.check_fg("Tertiary alcohol", product_smiles)
                    ):
                        reaction_sequence.append(("reduction", depth))
                        print(f"Detected reduction: Ester to alcohol at depth {depth}")

                    elif any(
                        checker.check_fg("Aldehyde", r) for r in reactants_smiles
                    ) and checker.check_fg("Primary alcohol", product_smiles):
                        reaction_sequence.append(("reduction", depth))
                        print(f"Detected reduction: Aldehyde to primary alcohol at depth {depth}")

                    elif any(
                        checker.check_fg("Ketone", r) for r in reactants_smiles
                    ) and checker.check_fg("Secondary alcohol", product_smiles):
                        reaction_sequence.append(("reduction", depth))
                        print(f"Detected reduction: Ketone to secondary alcohol at depth {depth}")

                    elif any(
                        checker.check_fg("Nitro group", r) for r in reactants_smiles
                    ) and checker.check_fg("Primary amine", product_smiles):
                        reaction_sequence.append(("reduction", depth))
                        print(f"Detected reduction: Nitro to amine at depth {depth}")

                    # Additional oxidation checks based on functional groups
                    elif any(
                        checker.check_fg("Primary alcohol", r) for r in reactants_smiles
                    ) and checker.check_fg("Aldehyde", product_smiles):
                        reaction_sequence.append(("oxidation", depth))
                        print(f"Detected oxidation: Primary alcohol to aldehyde at depth {depth}")

                    elif any(
                        checker.check_fg("Secondary alcohol", r) for r in reactants_smiles
                    ) and checker.check_fg("Ketone", product_smiles):
                        reaction_sequence.append(("oxidation", depth))
                        print(f"Detected oxidation: Secondary alcohol to ketone at depth {depth}")

                    elif any(
                        checker.check_fg("Primary alcohol", r) for r in reactants_smiles
                    ) and checker.check_fg("Carboxylic acid", product_smiles):
                        reaction_sequence.append(("oxidation", depth))
                        print(
                            f"Detected oxidation: Primary alcohol to carboxylic acid at depth {depth}"
                        )

                    elif any(
                        checker.check_fg("Aldehyde", r) for r in reactants_smiles
                    ) and checker.check_fg("Carboxylic acid", product_smiles):
                        reaction_sequence.append(("oxidation", depth))
                        print(f"Detected oxidation: Aldehyde to carboxylic acid at depth {depth}")

                    elif any(checker.check_fg("Primary amine", r) for r in reactants_smiles) and (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    ):
                        reaction_sequence.append(("oxidation", depth))
                        print(f"Detected oxidation: Amine to amide at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        if "children" in node:
            for child in node["children"]:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Count reduction and oxidation steps
    reduction_steps = sum(1 for step_type, _ in reaction_sequence if step_type == "reduction")
    oxidation_steps = sum(1 for step_type, _ in reaction_sequence if step_type == "oxidation")

    print(f"Reduction steps: {reduction_steps}")
    print(f"Oxidation steps: {oxidation_steps}")

    # Check if we have a sequence of reduction and oxidation steps
    has_sequence = False

    if len(reaction_sequence) >= 2:
        # Sort by depth (retrosynthetic order)
        reaction_sequence.sort(key=lambda x: x[1])

        # Check for sequences
        for i in range(len(reaction_sequence) - 1):
            if (
                reaction_sequence[i][0] == "reduction"
                and reaction_sequence[i + 1][0] == "oxidation"
            ) or (
                reaction_sequence[i][0] == "oxidation"
                and reaction_sequence[i + 1][0] == "reduction"
            ):
                has_sequence = True
                print(
                    f"Found reduction-oxidation sequence: {reaction_sequence[i][0]} at depth {reaction_sequence[i][1]} followed by {reaction_sequence[i+1][0]} at depth {reaction_sequence[i+1][1]}"
                )
                break

    return has_sequence and reduction_steps >= 1 and oxidation_steps >= 1
