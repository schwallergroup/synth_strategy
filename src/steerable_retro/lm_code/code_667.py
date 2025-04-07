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
    Detects a synthetic strategy where an azide intermediate is used to convert
    an alkyl halide to an amine, which is then protected with a Boc group.

    In retrosynthetic analysis:
    Boc-protected amine → Primary amine → Azide → Alkyl halide
    """
    # Initialize tracking variables
    steps_found = {"boc_protection": False, "azide_reduction": False, "azide_formation": False}

    # Track the depth of each step to ensure correct sequence
    step_depths = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # In retrosynthetic analysis, the product is the starting point
                # and reactants are what it's broken down into

                # Check for Boc protection (first step in retrosynthesis)
                boc_protection_reactions = [
                    "Boc amine protection",
                    "Boc amine protection explicit",
                    "Boc amine protection with Boc anhydride",
                    "Boc amine protection (ethyl Boc)",
                    "Boc amine protection of primary amine",
                ]

                if any(checker.check_reaction(rxn, rsmi) for rxn in boc_protection_reactions):
                    steps_found["boc_protection"] = True
                    step_depths["boc_protection"] = depth
                    print(f"Found Boc protection step at depth {depth}")
                elif checker.check_fg("Carbamic ester", product) and any(
                    checker.check_fg("Primary amine", reactant) for reactant in reactants
                ):
                    steps_found["boc_protection"] = True
                    step_depths["boc_protection"] = depth
                    print(
                        f"Found Boc protection step at depth {depth} (detected by functional groups)"
                    )

                # Check for azide reduction to amine
                if checker.check_reaction("Azide to amine reduction (Staudinger)", rsmi):
                    steps_found["azide_reduction"] = True
                    step_depths["azide_reduction"] = depth
                    print(f"Found azide reduction step at depth {depth}")
                elif any(
                    checker.check_fg("Azide", reactant) for reactant in reactants
                ) and checker.check_fg("Primary amine", product):
                    steps_found["azide_reduction"] = True
                    step_depths["azide_reduction"] = depth
                    print(
                        f"Found azide reduction step at depth {depth} (detected by functional groups)"
                    )

                # Check for azide formation from alkyl halide
                if checker.check_reaction("Formation of Azides from halogens", rsmi):
                    steps_found["azide_formation"] = True
                    step_depths["azide_formation"] = depth
                    print(f"Found azide formation step at depth {depth}")
                elif any(
                    checker.check_fg("Primary halide", reactant)
                    or checker.check_fg("Secondary halide", reactant)
                    or checker.check_fg("Tertiary halide", reactant)
                    for reactant in reactants
                ) and checker.check_fg("Azide", product):
                    steps_found["azide_formation"] = True
                    step_depths["azide_formation"] = depth
                    print(
                        f"Found azide formation step at depth {depth} (detected by functional groups)"
                    )

        # Traverse children (moving deeper in retrosynthesis)
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all steps were found
    all_steps_found = all(steps_found.values())

    # Check if steps are in the correct sequence in retrosynthesis:
    # Boc-protected amine → Primary amine → Azide → Alkyl halide
    correct_sequence = False
    if all_steps_found:
        if (
            step_depths["boc_protection"]
            < step_depths["azide_reduction"]
            < step_depths["azide_formation"]
        ):
            correct_sequence = True
            print("All steps found in correct sequence")
        else:
            print("Steps found but in incorrect sequence")

    return all_steps_found and correct_sequence
