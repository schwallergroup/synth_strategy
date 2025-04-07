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
    Detects if the synthetic route involves reductive amination (aldehyde/ketone + amine → amine).
    """
    reductive_amination_detected = False

    def dfs_traverse(node):
        nonlocal reductive_amination_detected

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a reductive amination reaction
            if (
                checker.check_reaction("reductive amination with aldehyde", rsmi)
                or checker.check_reaction("reductive amination with ketone", rsmi)
                or checker.check_reaction("reductive amination with alcohol", rsmi)
            ):
                print(f"Detected reductive amination: {rsmi}")
                reductive_amination_detected = True

            # If direct reaction check fails, check for components
            if not reductive_amination_detected:
                reactants_str = rsmi.split(">")[0]
                product_str = rsmi.split(">")[-1]
                reactants = reactants_str.split(".")

                # Check for aldehyde/ketone and amine in reactants
                has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
                has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
                has_amine = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                )

                # Check if product has amine but not imine (completed reduction)
                if (
                    has_amine
                    and (has_aldehyde or has_ketone)
                    and checker.check_fg("Secondary amine", product_str)
                    or checker.check_fg("Tertiary amine", product_str)
                ):
                    print(f"Detected possible reductive amination components: {rsmi}")
                    reductive_amination_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return reductive_amination_detected
