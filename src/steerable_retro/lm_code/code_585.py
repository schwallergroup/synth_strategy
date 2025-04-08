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
    Detects a strategy involving reduction of aromatic nitro to amine,
    followed by further functionalization of the amine.
    """
    # Track molecules and reactions in the synthesis path
    synthesis_paths = []

    def dfs_traverse(node, current_path=None):
        if current_path is None:
            current_path = []

        # Add current node to path
        current_path.append(node)

        # If this is a leaf node (starting material), save the path
        if node["type"] == "mol" and node.get("in_stock", False):
            synthesis_paths.append(current_path.copy())

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_path.copy())

    # Start traversal to collect all synthesis paths
    dfs_traverse(route)

    print(f"Found {len(synthesis_paths)} synthesis paths")

    # Check each synthesis path for the nitro reduction -> amine functionalization pattern
    for path_idx, path in enumerate(synthesis_paths):
        print(f"Analyzing path {path_idx+1}/{len(synthesis_paths)}")

        # Find amine functionalization and subsequent nitro reduction
        # Note: In retrosynthetic analysis, we see amine functionalization first, then nitro reduction
        amine_functionalization_indices = []
        nitro_reduction_indices = []

        for i in range(len(path) - 1):
            if path[i]["type"] == "reaction" and path[i + 1]["type"] == "mol":
                reaction_node = path[i]
                product_node = path[i + 1]

                # Skip if no reaction SMILES
                if "rsmi" not in reaction_node["metadata"]:
                    print(f"No rsmi in reaction at index {i}")
                    continue

                rsmi = reaction_node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Analyzing reaction at index {i}: {rsmi}")

                # Check for amine functionalization
                amine_functionalization_reactions = [
                    "N-alkylation of primary amines with alkyl halides",
                    "N-alkylation of secondary amines with alkyl halides",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine",
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine",
                    "Reductive amination with aldehyde",
                    "Reductive amination with ketone",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine",
                    "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine",
                    "Urea synthesis via isocyanate and primary amine",
                    "Urea synthesis via isocyanate and secondary amine",
                ]

                # Check for specific reductive amination or methylation patterns
                if checker.check_fg("Tertiary amine", product) and any(
                    checker.check_fg("Aniline", r)
                    or checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    for r in reactants
                ):
                    print(f"Found amine functionalization (methylation) at index {i}: {rsmi}")
                    amine_functionalization_indices.append(i)
                elif any(
                    checker.check_reaction(rxn, rsmi) for rxn in amine_functionalization_reactions
                ):
                    # Check if reactant has primary or secondary amine
                    if any(
                        checker.check_fg("Aniline", r)
                        or checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    ):
                        print(f"Found amine functionalization at index {i}: {rsmi}")
                        amine_functionalization_indices.append(i)

                # Check for nitro reduction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction at index {i}: {rsmi}")
                    nitro_reduction_indices.append(i)
                # Alternative check for nitro reduction
                elif any(checker.check_fg("Nitro group", r) for r in reactants) and (
                    checker.check_fg("Aniline", product)
                    or checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                ):
                    print(f"Found nitro reduction (alternative check) at index {i}: {rsmi}")
                    nitro_reduction_indices.append(i)

        # Check if we have both reactions in the correct order in the retrosynthetic direction
        # This means amine functionalization should come before nitro reduction
        for amine_idx in amine_functionalization_indices:
            for nitro_idx in nitro_reduction_indices:
                if (
                    nitro_idx > amine_idx
                ):  # In retrosynthetic direction, nitro reduction comes after amine functionalization
                    print(
                        f"Found complete pattern: amine functionalization at {amine_idx}, nitro reduction at {nitro_idx}"
                    )
                    return True

    print("No nitro reduction -> amine functionalization pattern found")
    return False
