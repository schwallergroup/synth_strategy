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
    Detects a synthetic strategy where a nitroaldol condensation is used to form a nitrostyrene,
    which is then reduced to a phenethylamine and further functionalized.
    """
    # Track if we found the key reactions
    found_nitroaldol = False
    found_nitroalkene_reduction = False
    found_nitro_reduction = False
    found_n_functionalization = False

    def dfs_traverse(node, depth=0):
        nonlocal found_nitroaldol, found_nitroalkene_reduction, found_nitro_reduction, found_n_functionalization

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitroaldol condensation (Henry reaction)
            if checker.check_reaction("Henry Reaction", rsmi):
                found_nitroaldol = True
                print(f"Found nitroaldol condensation at depth {depth}: {rsmi}")
            # Alternative check for nitroaldol if reaction check fails
            elif (
                any(checker.check_fg("Aldehyde", r) for r in reactants)
                and any(checker.check_fg("Nitro group", r) for r in reactants)
                and checker.check_fg("Nitro group", product)
                and "C=C" in product
            ):
                found_nitroaldol = True
                print(f"Found nitroaldol condensation (by FG) at depth {depth}: {rsmi}")

            # Check for nitroalkene reduction (either direct or stepwise)
            # First check for direct reduction of nitroalkene to phenethylamine
            if any(checker.check_fg("Nitro group", r) and "C=C" in r for r in reactants):
                if checker.check_fg("Primary amine", product) and not checker.check_fg(
                    "Nitro group", product
                ):
                    found_nitroalkene_reduction = True
                    print(f"Found nitroalkene reduction at depth {depth}: {rsmi}")

            # Check for nitro reduction to amine (part of stepwise reduction)
            if any(checker.check_fg("Nitro group", r) for r in reactants) and not any(
                checker.check_fg("Nitro group", r) for r in [product]
            ):
                if checker.check_fg("Primary amine", product):
                    found_nitro_reduction = True
                    print(f"Found nitro reduction at depth {depth}: {rsmi}")

                    # Try to check if this is specifically the "Reduction of nitro groups to amines" reaction
                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        print(f"Confirmed nitro reduction reaction type at depth {depth}")

            # Check for N-functionalization of the amine
            if any(checker.check_fg("Primary amine", r) for r in reactants):
                # Check for various N-functionalization products
                if (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                    or checker.check_fg("Urea", product)
                    or checker.check_fg("Sulfonamide", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                ):
                    found_n_functionalization = True
                    print(f"Found N-functionalization at depth {depth}: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Print summary of findings
    print(f"Nitroaldol condensation found: {found_nitroaldol}")
    print(f"Nitroalkene reduction found: {found_nitroalkene_reduction}")
    print(f"Nitro reduction found: {found_nitro_reduction}")
    print(f"N-functionalization found: {found_n_functionalization}")

    # Return True if we found the key components of this strategy
    # Either nitroaldol OR nitro_reduction is required, along with N-functionalization
    return (found_nitroaldol or found_nitro_reduction) and found_n_functionalization
