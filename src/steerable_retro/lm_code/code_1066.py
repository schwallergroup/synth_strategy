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
    This function detects a synthetic strategy with late-stage cyanation (final step)
    following multiple halogenation steps, starting from nitro reduction.
    """
    has_cyanation = False
    has_bromination = False
    has_nitro_reduction = False
    cyanation_depth = None
    bromination_depths = []
    nitro_reduction_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal has_cyanation, has_bromination, has_nitro_reduction
        nonlocal cyanation_depth, bromination_depths, nitro_reduction_depth

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for cyanation (halide to nitrile)
            if any(checker.check_fg("Nitrile", product) for _ in [1]) and any(
                checker.check_fg("Aromatic halide", r) for r in reactants
            ):
                has_cyanation = True
                cyanation_depth = depth
                print(f"Found cyanation at depth {depth}: {rsmi}")

            # Check for bromination reactions
            bromination_reactions = [
                "Aromatic bromination",
                "Bromination",
                "Wohl-Ziegler bromination benzyl primary",
                "Wohl-Ziegler bromination benzyl secondary",
                "Wohl-Ziegler bromination benzyl tertiary",
                "Wohl-Ziegler bromination allyl primary",
                "Wohl-Ziegler bromination allyl secondary",
                "Wohl-Ziegler bromination allyl tertiary",
            ]

            if any(checker.check_reaction(rxn, rsmi) for rxn in bromination_reactions) or (
                any(checker.check_fg("Aromatic halide", product) for _ in [1])
                and not any(checker.check_fg("Aromatic halide", r) for r in reactants)
            ):
                has_bromination = True
                bromination_depths.append(depth)
                print(f"Found bromination at depth {depth}: {rsmi}")

            # Check for nitro reduction (NO2 to NH2)
            if (
                checker.check_reaction("Reduction of nitro groups to amines", rsmi)
                or (
                    any(checker.check_fg("Aniline", product) for _ in [1])
                    and any(checker.check_fg("Nitro group", r) for r in reactants)
                )
                or (
                    any(checker.check_fg("Primary amine", product) for _ in [1])
                    and any(checker.check_fg("Nitro group", r) for r in reactants)
                )
            ):
                has_nitro_reduction = True
                nitro_reduction_depth = depth
                print(f"Found nitro reduction at depth {depth}: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if the strategy criteria are met
    strategy_present = (
        has_cyanation
        and has_bromination
        and len(bromination_depths) >= 1
        and cyanation_depth <= 1  # Cyanation is a late-stage step (depth 0 or 1)
    )

    # If nitro reduction is present, ensure it occurs earlier than cyanation
    if has_nitro_reduction:
        strategy_present = strategy_present and (nitro_reduction_depth > cyanation_depth)

    print(f"Strategy detection result: {strategy_present}")
    print(f"Cyanation depth: {cyanation_depth}")
    print(f"Bromination depths: {bromination_depths}")
    print(f"Nitro reduction depth: {nitro_reduction_depth}")

    return strategy_present
