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
    This function detects if the synthesis begins with a nitroaldol condensation (C-C bond formation).
    Early stage corresponds to high depth in the retrosynthetic tree.
    """
    # Track if we found a Henry reaction and at what depth
    found_henry_reaction = False
    henry_reaction_depth = 0
    max_depth_found = 0

    def is_nitroaldol_condensation(rsmi):
        """Helper function to identify nitroaldol condensations based on functional groups"""
        try:
            # Split reaction SMILES into reactants and products
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if any reactant has a nitro group
            has_nitro_reactant = any(checker.check_fg("Nitro group", r) for r in reactants_smiles)

            # Check if any reactant has an aldehyde
            has_aldehyde_reactant = any(checker.check_fg("Aldehyde", r) for r in reactants_smiles)

            # If both nitro and aldehyde are present in reactants, it's likely a nitroaldol condensation
            if has_nitro_reactant and has_aldehyde_reactant:
                print(f"Detected nitroaldol condensation pattern in: {rsmi}")
                return True

            return False
        except Exception as e:
            print(f"Error in is_nitroaldol_condensation: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal found_henry_reaction, henry_reaction_depth, max_depth_found

        # Update max depth
        max_depth_found = max(max_depth_found, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            # Check if this is a Henry Reaction using the checker or by functional group pattern
            if checker.check_reaction("Henry Reaction", rsmi) or is_nitroaldol_condensation(rsmi):
                print(f"Found nitroaldol condensation at depth {depth}")
                found_henry_reaction = True
                henry_reaction_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Traverse the route to find Henry reactions and max depth
    dfs_traverse(route)

    # Determine if the Henry reaction is in the early stage (high depth)
    # Consider it early stage if it's within 3 levels of the maximum depth
    early_stage_henry = found_henry_reaction and henry_reaction_depth >= max_depth_found - 3

    print(f"Max depth found: {max_depth_found}")
    print(
        f"Nitroaldol condensation depth: {henry_reaction_depth if found_henry_reaction else 'Not found'}"
    )
    print(f"Early nitroaldol condensation strategy: {early_stage_henry}")

    return early_stage_henry
