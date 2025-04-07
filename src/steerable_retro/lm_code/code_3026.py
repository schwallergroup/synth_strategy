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
    Detects if the synthesis uses multiple Friedel-Crafts reactions in sequence.
    """
    friedel_crafts_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check if this is a Friedel-Crafts reaction using the checker function
            is_friedel_crafts_acylation = checker.check_reaction("Friedel-Crafts acylation", rsmi)
            is_friedel_crafts_alkylation = checker.check_reaction("Friedel-Crafts alkylation", rsmi)
            is_friedel_crafts_alkylation_halide = checker.check_reaction(
                "Friedel-Crafts alkylation with halide", rsmi
            )

            print(f"Is Friedel-Crafts acylation: {is_friedel_crafts_acylation}")
            print(f"Is Friedel-Crafts alkylation: {is_friedel_crafts_alkylation}")
            print(
                f"Is Friedel-Crafts alkylation with halide: {is_friedel_crafts_alkylation_halide}"
            )

            # If direct check fails, try pattern-based detection
            is_friedel_crafts = (
                is_friedel_crafts_acylation
                or is_friedel_crafts_alkylation
                or is_friedel_crafts_alkylation_halide
            )

            if not is_friedel_crafts:
                # Parse reactants and products
                try:
                    reactants_part = rsmi.split(">")[0]
                    reagents_part = rsmi.split(">")[1] if len(rsmi.split(">")) > 2 else ""
                    product_part = rsmi.split(">")[-1]

                    reactants = reactants_part.split(".")
                    reagents = reagents_part.split(".") if reagents_part else []

                    # Check for Lewis acid catalysts (common in Friedel-Crafts)
                    lewis_acid_present = any(
                        re.search(r"Al|FeCl|AlCl|BF3|TiCl", r) for r in reagents
                    )

                    # Check for acyl halides (for acylation)
                    acyl_halide_present = any(checker.check_fg("Acyl halide", r) for r in reactants)

                    # Check for aromatic compounds
                    aromatic_present = any("c" in r.lower() for r in reactants)

                    # Check for alkyl halides (for alkylation)
                    alkyl_halide_present = any(
                        checker.check_fg("Primary halide", r)
                        or checker.check_fg("Secondary halide", r)
                        or checker.check_fg("Tertiary halide", r)
                        for r in reactants
                    )

                    # Detect Friedel-Crafts acylation
                    if (
                        acyl_halide_present
                        and aromatic_present
                        and (
                            lewis_acid_present
                            or "AlCl3" in reagents_part
                            or "[Al+3]" in reagents_part
                        )
                    ):
                        is_friedel_crafts = True
                        print(f"Detected Friedel-Crafts acylation pattern at depth {depth}")

                    # Detect Friedel-Crafts alkylation
                    elif (
                        alkyl_halide_present
                        and aromatic_present
                        and (
                            lewis_acid_present
                            or "AlCl3" in reagents_part
                            or "[Al+3]" in reagents_part
                        )
                    ):
                        is_friedel_crafts = True
                        print(f"Detected Friedel-Crafts alkylation pattern at depth {depth}")

                except Exception as e:
                    print(f"Error analyzing reaction: {e}")

            if is_friedel_crafts:
                friedel_crafts_reactions.append(depth)
                print(f"Found Friedel-Crafts reaction at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"All Friedel-Crafts reactions found at depths: {friedel_crafts_reactions}")

    # Check if we have at least 2 Friedel-Crafts reactions
    if len(friedel_crafts_reactions) >= 2:
        # Check if they are in sequence (adjacent depths or with only 1 step between)
        friedel_crafts_reactions.sort()
        for i in range(len(friedel_crafts_reactions) - 1):
            if friedel_crafts_reactions[i + 1] - friedel_crafts_reactions[i] <= 2:
                print(
                    f"Found Friedel-Crafts sequence at depths {friedel_crafts_reactions[i]} and {friedel_crafts_reactions[i+1]}"
                )
                return True

    return False
