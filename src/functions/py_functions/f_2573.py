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
    Detects a synthetic strategy involving late-stage C-C bond formation via coupling
    followed by sequential oxidation steps (vinyl → aldehyde → carboxylic acid)
    and ending with functional group activation (acid → acid chloride).
    """
    # Initialize tracking variables
    reactions_with_depth = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_part = rsmi.split(">")[0]
                products_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")
                product = products_part

                # Check for C-C coupling reactions (Heck, Suzuki, Stille, etc.)
                cc_coupling = False
                if (
                    checker.check_reaction("Heck terminal vinyl", rsmi)
                    or checker.check_reaction("Heck_terminal_vinyl", rsmi)
                    or checker.check_reaction(
                        "Suzuki coupling with boronic acids", rsmi
                    )
                    or checker.check_reaction(
                        "Suzuki coupling with boronic esters", rsmi
                    )
                    or checker.check_reaction("Suzuki", rsmi)
                    or checker.check_reaction("Stille reaction_vinyl", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Stille", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Negishi", rsmi)
                    or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                    or checker.check_reaction("Sonogashira acetylene_aryl halide", rsmi)
                ):

                    # Verify a new carbon-carbon bond is formed
                    if (
                        checker.check_fg("Vinyl", product)
                        or checker.check_fg("Alkyne", product)
                        or checker.check_fg("Allyl", product)
                    ):
                        cc_coupling = True
                        print(f"Detected C-C coupling at depth {depth}")
                        reactions_with_depth.append(("cc_coupling", depth))

                # Check for vinyl to aldehyde oxidation
                if (
                    checker.check_reaction("Alkene oxidation to aldehyde", rsmi)
                    or checker.check_reaction(
                        "Oxidation of alkene to carboxylic acid", rsmi
                    )
                    or (
                        any(checker.check_fg("Vinyl", r) for r in reactants)
                        and checker.check_fg("Aldehyde", product)
                        and not any(checker.check_fg("Aldehyde", r) for r in reactants)
                    )
                ):
                    print(f"Detected vinyl to aldehyde oxidation at depth {depth}")
                    reactions_with_depth.append(("vinyl_to_aldehyde", depth))

                # Check for aldehyde to carboxylic acid oxidation
                if checker.check_reaction(
                    "Oxidation of aldehydes to carboxylic acids", rsmi
                ) or (
                    any(checker.check_fg("Aldehyde", r) for r in reactants)
                    and checker.check_fg("Carboxylic acid", product)
                    and not any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                ):
                    print(
                        f"Detected aldehyde to carboxylic acid oxidation at depth {depth}"
                    )
                    reactions_with_depth.append(("aldehyde_to_acid", depth))

                # Check for acid to acid chloride activation
                if (
                    checker.check_reaction("Acyl chlorides from alcohols", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with ammonia to amide", rsmi
                    )
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        rsmi,
                    )
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or (
                        any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                        and checker.check_fg("Acyl halide", product)
                        and not any(
                            checker.check_fg("Acyl halide", r) for r in reactants
                        )
                    )
                ):
                    print(
                        f"Detected carboxylic acid to acid chloride activation at depth {depth}"
                    )
                    reactions_with_depth.append(("acid_to_acid_chloride", depth))

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (ascending - lower depth is later in synthesis)
    reactions_with_depth.sort(key=lambda x: x[1])

    # Extract just the reaction types in order
    reaction_sequence = [r[0] for r in reactions_with_depth]
    print(f"Reaction sequence (retrosynthetic): {reaction_sequence}")

    # Check if we have the required reactions
    has_cc_coupling = "cc_coupling" in reaction_sequence
    has_vinyl_to_aldehyde = "vinyl_to_aldehyde" in reaction_sequence
    has_aldehyde_to_acid = "aldehyde_to_acid" in reaction_sequence
    has_acid_to_acid_chloride = "acid_to_acid_chloride" in reaction_sequence

    # Check if we have at least two oxidation steps
    has_sequential_oxidation = has_vinyl_to_aldehyde and has_aldehyde_to_acid

    # In retrosynthetic analysis, the sequence should be:
    # acid_to_acid_chloride → aldehyde_to_acid → vinyl_to_aldehyde → cc_coupling
    # (which is the reverse of the forward synthesis)

    correct_sequence = False

    # We need at least the sequential oxidation steps
    if has_sequential_oxidation:
        try:
            aldehyde_to_acid_idx = reaction_sequence.index("aldehyde_to_acid")
            vinyl_to_aldehyde_idx = reaction_sequence.index("vinyl_to_aldehyde")

            # Check if oxidation steps are in correct order
            if aldehyde_to_acid_idx < vinyl_to_aldehyde_idx:
                # If we have C-C coupling, it should come after vinyl_to_aldehyde
                if not has_cc_coupling or (
                    has_cc_coupling
                    and vinyl_to_aldehyde_idx < reaction_sequence.index("cc_coupling")
                ):
                    correct_sequence = True
        except ValueError:
            pass

    # For the full strategy, we need C-C coupling and sequential oxidation in correct order
    # But we'll also accept just the sequential oxidation as a partial match
    strategy_present = has_sequential_oxidation and correct_sequence

    if strategy_present:
        if has_cc_coupling:
            print(
                "Detected complete sequential oxidation strategy with late-stage C-C bond formation"
            )
        else:
            print("Detected sequential oxidation strategy (without C-C coupling)")
    else:
        print("Strategy not detected")
        if not has_sequential_oxidation:
            print("Missing sequential oxidation steps")
        if not correct_sequence:
            print("Incorrect reaction sequence")

    return strategy_present
