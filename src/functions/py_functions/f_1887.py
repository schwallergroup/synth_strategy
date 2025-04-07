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
    Detects amide formation as part of fragment preparation before key coupling steps.
    """
    # Track transformations and their depths
    amide_formation = False
    coupling_reaction = False

    amide_depth = -1
    coupling_depth = -1

    # Lists of relevant reaction types
    amide_formation_reactions = [
        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
        "Carboxylic acid with primary amine to amide",
        "Ester with primary amine to amide",
        "Ester with ammonia to amide",
        "Ester with secondary amine to amide",
        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
        "Acyl chloride with secondary amine to amide",
        "Acyl chloride with ammonia to amide",
        "Schotten-Baumann_amide",
    ]

    coupling_reactions = [
        "Suzuki coupling with boronic acids",
        "Suzuki coupling with boronic esters",
        "Suzuki coupling with boronic acids OTf",
        "Suzuki coupling with boronic esters OTf",
        "Suzuki",
        "Negishi coupling",
        "Negishi",
        "Stille reaction_aryl",
        "Stille reaction_vinyl",
        "Stille",
        "Heck terminal vinyl",
        "Heck_terminal_vinyl",
        "Sonogashira alkyne_aryl halide",
        "Sonogashira acetylene_aryl halide",
        "Hiyama-Denmark Coupling",
        "Kumada cross-coupling",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation, coupling_reaction
        nonlocal amide_depth, coupling_depth

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check for amide formation using predefined reaction types
                amide_found = False
                for rxn_type in amide_formation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        amide_found = True
                        amide_formation = True
                        amide_depth = depth
                        print(f"Found amide formation ({rxn_type}) at depth {depth}")
                        break

                # If no predefined amide formation reaction was found, check manually
                if not amide_found:
                    # Extract reactants and product
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if product contains amide group
                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        ):
                            # Check if reactants contain carboxylic acid and amine
                            has_carboxylic_acid = any(
                                checker.check_fg("Carboxylic acid", r)
                                for r in reactants
                            )
                            has_amine = any(
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Secondary amine", r)
                                for r in reactants
                            )

                            if has_carboxylic_acid and has_amine:
                                amide_formation = True
                                amide_depth = depth
                                print(
                                    f"Found amide formation (manual check) at depth {depth}"
                                )

                        # Also check for acyl chloride + amine → amide
                        if not amide_formation and (
                            checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Primary amide", product)
                        ):
                            has_acyl_halide = any(
                                checker.check_fg("Acyl halide", r) for r in reactants
                            )
                            has_amine = any(
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Secondary amine", r)
                                for r in reactants
                            )

                            if has_acyl_halide and has_amine:
                                amide_formation = True
                                amide_depth = depth
                                print(
                                    f"Found amide formation (acyl halide + amine) at depth {depth}"
                                )
                    except Exception as e:
                        print(f"Error in manual amide check: {e}")

                # Check for coupling reaction using predefined reaction types
                coupling_found = False
                for rxn_type in coupling_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        coupling_found = True
                        coupling_reaction = True
                        coupling_depth = depth
                        print(f"Found coupling reaction ({rxn_type}) at depth {depth}")
                        break

                # If no predefined coupling reaction was found, check manually
                if not coupling_found:
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check for boronic acid/ester in reactants (Suzuki)
                        has_boronic = any(
                            checker.check_fg("Boronic acid", r)
                            or checker.check_fg("Boronic ester", r)
                            for r in reactants
                        )
                        has_halide = any(
                            checker.check_fg("Aromatic halide", r) for r in reactants
                        )

                        if has_boronic and has_halide:
                            coupling_reaction = True
                            coupling_depth = depth
                            print(
                                f"Found coupling reaction (manual Suzuki check) at depth {depth}"
                            )
                    except Exception as e:
                        print(f"Error in manual coupling check: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if amide formation occurs before coupling (higher depth = earlier in synthesis)
    amide_before_coupling = (
        amide_formation and coupling_reaction and amide_depth > coupling_depth
    )

    print(f"Amide formation: {amide_formation} at depth {amide_depth}")
    print(f"Coupling reaction: {coupling_reaction} at depth {coupling_depth}")
    print(f"Amide formation in fragment preparation detected: {amide_before_coupling}")

    return amide_before_coupling
