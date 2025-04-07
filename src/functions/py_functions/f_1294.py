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
    This function detects reductive amination (aldehyde + amine → secondary amine).
    """
    reductive_amination_detected = False

    def dfs_traverse(node):
        nonlocal reductive_amination_detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if the reaction is directly identified as reductive amination
            is_reductive_amination = checker.check_reaction(
                "Reductive amination with aldehyde", rsmi
            )

            # Alternatively, check for the presence of required functional groups
            has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
            has_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
            forms_secondary_amine = checker.check_fg(
                "Secondary amine", product
            ) and not checker.check_fg("Primary amine", product)

            # Check for reductive amination with ketone as well
            is_reductive_amination_ketone = checker.check_reaction(
                "Reductive amination with ketone", rsmi
            )
            has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)

            if (
                is_reductive_amination
                or is_reductive_amination_ketone
                or (
                    (has_aldehyde or has_ketone) and has_amine and forms_secondary_amine
                )
            ):
                reductive_amination_detected = True
                print(f"Detected reductive amination: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return reductive_amination_detected
