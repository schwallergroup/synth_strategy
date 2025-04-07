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
    Detects a nitro reduction to amine followed by sulfonamide formation sequence.

    In retrosynthetic traversal, we'll encounter sulfonamide formation first,
    then nitro reduction.
    """
    # Track if we found the sequence
    sequence_found = False

    # Track the reactions in sequence (in retrosynthetic order)
    sulfonamide_formation_found = False
    nitro_reduction_found = False

    # Track the molecule containing the amine that links the two reactions
    amine_containing_mol = None

    def dfs_traverse(node, depth=0):
        nonlocal sequence_found, sulfonamide_formation_found, nitro_reduction_found, amine_containing_mol

        if sequence_found:
            return  # Stop traversal if sequence already found

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Depth {depth}, Examining reaction: {rsmi}")

                # In retrosynthetic traversal, we encounter sulfonamide formation first
                if not sulfonamide_formation_found:
                    # Check if this is a sulfonamide formation reaction
                    is_sulfonamide_reaction = checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                    ) or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                    )

                    # Alternative check: look for amine reactant and sulfonamide product
                    has_amine_reactant = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )
                    has_sulfonamide_product = checker.check_fg("Sulfonamide", product)

                    if is_sulfonamide_reaction or (has_amine_reactant and has_sulfonamide_product):
                        print(f"Found sulfonamide formation at depth {depth}")
                        sulfonamide_formation_found = True

                        # Find which reactant contains the amine
                        for reactant in reactants:
                            if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                                "Secondary amine", reactant
                            ):
                                amine_containing_mol = reactant
                                print(f"Amine-containing molecule: {amine_containing_mol}")
                                break

                # After finding sulfonamide formation, look for nitro reduction
                elif (
                    sulfonamide_formation_found
                    and not nitro_reduction_found
                    and amine_containing_mol is not None
                ):
                    # Check if this is a nitro reduction reaction
                    is_nitro_reduction = checker.check_reaction(
                        "Reduction of nitro groups to amines", rsmi
                    )

                    # Alternative check: look for nitro reactant and amine product
                    has_nitro_reactant = any(checker.check_fg("Nitro group", r) for r in reactants)
                    has_amine_product = checker.check_fg(
                        "Primary amine", product
                    ) or checker.check_fg("Secondary amine", product)

                    if is_nitro_reduction or (has_nitro_reactant and has_amine_product):
                        print(f"Found potential nitro reduction at depth {depth}")

                        # Check if the product of this reaction is the same as the amine used in sulfonamide formation
                        # This can be a direct string comparison or a more flexible check
                        if product == amine_containing_mol or (
                            checker.check_fg("Primary amine", product)
                            and checker.check_fg("Primary amine", amine_containing_mol)
                        ):
                            print("Confirmed nitro reduction to the same amine")
                            nitro_reduction_found = True
                            sequence_found = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Sequence found: {sequence_found}")
    print(f"Sulfonamide formation found: {sulfonamide_formation_found}")
    print(f"Nitro reduction found: {nitro_reduction_found}")

    return sequence_found
