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
    This function detects a synthetic strategy involving the conversion of an alcohol
    to an amine through a sequence of functional group interconversions:
    alcohol → mesylate/triflate → azide → amine
    """
    # Track transformations and their sequence
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for alcohol to mesylate/triflate transformation
                if any(
                    checker.check_fg("Primary alcohol", r)
                    or checker.check_fg("Secondary alcohol", r)
                    or checker.check_fg("Tertiary alcohol", r)
                    for r in reactants
                ) and (
                    checker.check_fg("Mesylate", product) or checker.check_fg("Triflate", product)
                ):
                    print(f"Found alcohol to mesylate/triflate transformation at depth {depth}")
                    transformations.append(("alcohol_to_leaving_group", depth))

                # Check for mesylate/triflate to azide transformation
                if any(
                    checker.check_fg("Mesylate", r) or checker.check_fg("Triflate", r)
                    for r in reactants
                ) and checker.check_fg("Azide", product):
                    print(f"Found mesylate/triflate to azide transformation at depth {depth}")
                    transformations.append(("leaving_group_to_azide", depth))

                # Check for azide to amine transformation
                if any(checker.check_fg("Azide", r) for r in reactants) and checker.check_fg(
                    "Primary amine", product
                ):
                    print(f"Found azide to amine transformation at depth {depth}")
                    transformations.append(("azide_to_amine", depth))

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have all three transformations
    transformation_types = [t[0] for t in transformations]

    if not all(
        t in transformation_types
        for t in ["alcohol_to_leaving_group", "leaving_group_to_azide", "azide_to_amine"]
    ):
        print("Not all required transformations found")
        return False

    # Sort transformations by depth to check sequence
    sorted_transformations = sorted(transformations, key=lambda x: x[1])
    sorted_types = [t[0] for t in sorted_transformations]

    # Check if the transformations appear in the correct sequence in the route
    # Note: In retrosynthesis, the sequence is reversed (higher depth = earlier in synthesis)
    correct_sequence = False

    for i in range(len(sorted_types) - 2):
        if (
            sorted_types[i] == "azide_to_amine"
            and sorted_types[i + 1] == "leaving_group_to_azide"
            and sorted_types[i + 2] == "alcohol_to_leaving_group"
        ):
            correct_sequence = True
            break

    print(f"Transformations in sequence: {correct_sequence}")
    return correct_sequence
