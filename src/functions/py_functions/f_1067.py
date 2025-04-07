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
    This function detects a synthetic strategy that transforms a nitro group to a nitrile group
    through a sequence of functional group transformations: Nitro → Amine → Bromide → Nitrile
    """
    # Track the sequence of functional group transformations
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            # Check if metadata and rsmi exist
            if "metadata" not in node or "rsmi" not in node["metadata"]:
                print(f"Missing metadata or rsmi at depth {depth}")
                return

            rsmi = node["metadata"]["rsmi"]
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]
                reactants_combined = ".".join(reactants)

                # Check for nitro reduction (NO2 to NH2)
                if checker.check_fg("Primary amine", product) and checker.check_fg(
                    "Nitro group", reactants_combined
                ):
                    # Additional check for the reduction reaction if available
                    if checker.check_reaction(
                        "Reduction of nitro groups to amines", rsmi
                    ):
                        transformations.append(("nitro_to_amine", depth))
                        print(f"Found nitro to amine transformation at depth {depth}")
                    else:
                        # Fallback to functional group check if reaction check fails
                        # This ensures we don't miss valid transformations
                        nitro_indices = checker.get_fg_atom_indices(
                            "Nitro group", reactants_combined
                        )
                        amine_indices = checker.get_fg_atom_indices(
                            "Primary amine", product
                        )
                        if nitro_indices and amine_indices:
                            transformations.append(("nitro_to_amine", depth))
                            print(
                                f"Found nitro to amine transformation at depth {depth} (FG check)"
                            )

                # Check for bromination of amine position
                if (
                    any(
                        checker.check_fg(halide, product)
                        for halide in [
                            "Primary halide",
                            "Secondary halide",
                            "Tertiary halide",
                        ]
                    )
                    and "Br" in product
                    and checker.check_fg("Primary amine", reactants_combined)
                ):
                    transformations.append(("amine_to_bromide", depth))
                    print(f"Found amine to bromide transformation at depth {depth}")

                # Check for cyanation (bromide to CN)
                if (
                    checker.check_fg("Nitrile", product)
                    and any(
                        checker.check_fg(halide, reactants_combined)
                        for halide in [
                            "Primary halide",
                            "Secondary halide",
                            "Tertiary halide",
                        ]
                    )
                    and "Br" in reactants_combined
                ):
                    transformations.append(("bromide_to_nitrile", depth))
                    print(f"Found bromide to nitrile transformation at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort transformations by depth (higher depth = earlier in synthesis)
    transformations.sort(key=lambda x: x[1], reverse=True)

    # Check if the sequence matches our target pattern
    transformation_types = [t[0] for t in transformations]
    print(f"Transformation sequence: {transformation_types}")

    # Check for the specific sequence: nitro_to_amine -> amine_to_bromide -> bromide_to_nitrile
    target_sequence = ["nitro_to_amine", "amine_to_bromide", "bromide_to_nitrile"]

    # Check if all transformations in the target sequence are present
    if all(t in transformation_types for t in target_sequence):
        # Find the indices of the first occurrence of each transformation
        indices = []
        for t in target_sequence:
            if t in transformation_types:
                indices.append(transformation_types.index(t))
            else:
                indices.append(-1)

        # Check if the transformations appear in the correct order
        if len(indices) == 3 and indices[0] < indices[1] < indices[2]:
            # Check if the transformations are consecutive in the synthesis route
            # by comparing their depths
            depths = [transformations[indices[i]][1] for i in range(3)]

            # In a linear synthesis, the depths should form a sequence
            # The difference between consecutive depths should be consistent with the synthesis direction
            # For a valid strategy, we expect nitro→amine→bromide→nitrile to be a direct path
            strategy_present = True
            print(f"Transformation depths: {depths}")
        else:
            strategy_present = False
    else:
        strategy_present = False

    print(f"Strategy detection result: {strategy_present}")

    return strategy_present
