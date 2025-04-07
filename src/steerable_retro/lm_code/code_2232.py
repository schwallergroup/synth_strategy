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
    This function detects a linear synthesis route with a specific functional group
    progression: carboxylic acid → Weinreb amide → ketone → amide → amine → alcohol
    """
    # Track functional groups observed in order
    functional_group_sequence = []

    def dfs_traverse(node, depth=0):
        nonlocal functional_group_sequence

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants_part = rsmi.split(">")[0]
                    product = rsmi.split(">")[-1]

                    # Validate SMILES
                    reactants_mol = [Chem.MolFromSmiles(r) for r in reactants_part.split(".") if r]
                    product_mol = Chem.MolFromSmiles(product)

                    if not all(reactants_mol) or not product_mol:
                        print(f"Invalid SMILES at depth {depth}: {rsmi}")
                        return

                    reactants = reactants_part.split(".")

                    # Check for specific functional group transformations
                    # Carboxylic acid to Weinreb amide
                    if checker.check_fg("Carboxylic acid", reactants_part) and checker.check_fg(
                        "Tertiary amide", product
                    ):
                        functional_group_sequence.append(("carboxylic_acid_to_weinreb", depth))
                        print(f"Found carboxylic acid → Weinreb amide at depth {depth}")

                    # Weinreb amide to ketone
                    elif checker.check_reaction("Ketone from Weinreb amide", rsmi) or (
                        checker.check_fg("Tertiary amide", reactants_part)
                        and checker.check_fg("Ketone", product)
                    ):
                        functional_group_sequence.append(("weinreb_to_ketone", depth))
                        print(f"Found Weinreb amide → ketone at depth {depth}")

                    # Ketone to amide
                    elif checker.check_fg("Ketone", reactants_part) and (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    ):
                        functional_group_sequence.append(("ketone_to_amide", depth))
                        print(f"Found ketone → amide at depth {depth}")

                    # Amide to amine (reduction of amides to amines)
                    elif (
                        checker.check_reaction("Reduction of primary amides to amines", rsmi)
                        or checker.check_reaction("Reduction of secondary amides to amines", rsmi)
                        or checker.check_reaction("Reduction of tertiary amides to amines", rsmi)
                        or (
                            (
                                checker.check_fg("Primary amide", reactants_part)
                                or checker.check_fg("Secondary amide", reactants_part)
                                or checker.check_fg("Tertiary amide", reactants_part)
                            )
                            and checker.check_fg("Primary amine", product)
                        )
                    ):
                        functional_group_sequence.append(("amide_to_amine", depth))
                        print(f"Found amide → amine at depth {depth}")

                    # Amine to alcohol
                    elif checker.check_fg("Primary amine", reactants_part) and (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    ):
                        functional_group_sequence.append(("amine_to_alcohol", depth))
                        print(f"Found amine → alcohol at depth {depth}")

                    # Ketone to alcohol (alternative path)
                    elif (
                        checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                        or checker.check_reaction(
                            "Reduction of aldehydes and ketones to alcohols", rsmi
                        )
                        or (
                            checker.check_fg("Ketone", reactants_part)
                            and (
                                checker.check_fg("Secondary alcohol", product)
                                or checker.check_fg("Primary alcohol", product)
                            )
                        )
                    ):
                        functional_group_sequence.append(("ketone_to_alcohol", depth))
                        print(f"Found ketone → alcohol at depth {depth}")
                except Exception as e:
                    print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth (higher depth = earlier in synthesis)
    functional_group_sequence.sort(key=lambda x: x[1], reverse=True)

    # Extract just the transformation types
    transformations = [fg[0] for fg in functional_group_sequence]
    print(f"Functional group progression: {transformations}")

    # Define the expected progression patterns
    expected_progressions = [
        # Full progression
        [
            "carboxylic_acid_to_weinreb",
            "weinreb_to_ketone",
            "ketone_to_amide",
            "amide_to_amine",
            "amine_to_alcohol",
        ],
        # Alternative with direct ketone to alcohol
        ["carboxylic_acid_to_weinreb", "weinreb_to_ketone", "ketone_to_alcohol"],
        # Partial progressions (at least 2 steps)
        ["carboxylic_acid_to_weinreb", "weinreb_to_ketone"],
        ["weinreb_to_ketone", "ketone_to_amide"],
        ["ketone_to_amide", "amide_to_amine"],
        ["amide_to_amine", "amine_to_alcohol"],
        ["weinreb_to_ketone", "ketone_to_alcohol"],
    ]

    # Check if any of the expected progressions is a subsequence of our transformations
    def is_subsequence(subseq, seq):
        i, j = 0, 0
        while i < len(subseq) and j < len(seq):
            if subseq[i] == seq[j]:
                i += 1
            j += 1
        return i == len(subseq)

    has_progression = any(is_subsequence(prog, transformations) for prog in expected_progressions)

    # Require at least 2 transformations
    has_progression = has_progression and len(functional_group_sequence) >= 2

    print(f"Strategy detection result: {has_progression}")
    return has_progression
