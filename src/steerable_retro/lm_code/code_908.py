#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    Detects a synthetic pathway with multiple nitrogen functional group interconversions
    (amine → azide → amine → sulfonamide → amine).
    """
    # Track the sequence of transformations
    transformation_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Depth {depth}, Analyzing reaction: {rsmi}")

                # In retrosynthesis, product is the starting material and reactants are the precursors
                # So we need to check transformations in reverse

                # Check for azide to amine transformation (retrosynthetically: amine in product, azide in reactants)
                if (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                ) and any(checker.check_fg("Azide", r) for r in reactants):
                    transformation_sequence.append("azide_to_amine")
                    print(f"Depth {depth}, Detected azide to amine transformation")

                # Check for amine to azide transformation (retrosynthetically: azide in product, amine in reactants)
                if checker.check_fg("Azide", product) and any(
                    checker.check_fg("Primary amine", r) for r in reactants
                ):
                    transformation_sequence.append("amine_to_azide")
                    print(f"Depth {depth}, Detected amine to azide transformation")

                # Check for sulfonamide to amine transformation (retrosynthetically: amine in product, sulfonamide in reactants)
                if (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                ) and any(checker.check_fg("Sulfonamide", r) for r in reactants):
                    transformation_sequence.append("sulfonamide_to_amine")
                    print(f"Depth {depth}, Detected sulfonamide to amine transformation")

                # Check for amine to sulfonamide transformation (retrosynthetically: sulfonamide in product, amine in reactants)
                if checker.check_fg("Sulfonamide", product) and any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                ):
                    transformation_sequence.append("amine_to_sulfonamide")
                    print(f"Depth {depth}, Detected amine to sulfonamide transformation")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Transformation sequence: {transformation_sequence}")

    # Check for the required sequence of transformations
    # We're looking for amine → azide → amine → sulfonamide → amine
    # In retrosynthesis, this would be: amine_to_azide, azide_to_amine, amine_to_sulfonamide, sulfonamide_to_amine

    # First, check if we have at least 3 of the required transformations
    required_transformations = [
        "amine_to_azide",
        "azide_to_amine",
        "amine_to_sulfonamide",
        "sulfonamide_to_amine",
    ]

    # Count how many of the required transformations are found
    count = sum(1 for t in required_transformations if t in transformation_sequence)

    # We need at least 3 different types of transformations
    if count >= 3:
        print(f"Found {count} of the required transformations")

        # Now check for sequences of 3 consecutive transformations
        for i in range(len(transformation_sequence) - 2):
            # Check for amine → azide → amine sequence
            if transformation_sequence[i : i + 3] == [
                "amine_to_azide",
                "azide_to_amine",
                "amine_to_sulfonamide",
            ] or transformation_sequence[i : i + 3] == [
                "azide_to_amine",
                "amine_to_sulfonamide",
                "sulfonamide_to_amine",
            ]:
                print(f"Found required transformation sequence at position {i}")
                return True

            # Also check for partial sequences if we have all 4 transformations but not in perfect order
            if (
                count == 4
                and set(transformation_sequence[i : i + 3]).issubset(set(required_transformations))
                and len(set(transformation_sequence[i : i + 3])) == 3
            ):
                print(f"Found 3 different required transformations in sequence at position {i}")
                return True

    # If we have all 4 transformations but not in sequence, still return True
    if count == 4:
        print("Found all 4 required transformations, but not in perfect sequence")
        return True

    # If we have 3 transformations and they're all different, return True
    if count == 3 and len(set(transformation_sequence)) >= 3:
        print("Found 3 different required transformations")
        return True

    return False
