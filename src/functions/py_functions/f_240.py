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
    Detects if the synthesis follows the strategy of early Grignard reaction for tertiary alcohol
    formation followed by late-stage aryl halide to ester conversion.
    """
    print(f"Analyzing route for early Grignard, late ester strategy...")
    grignard_depths = []
    ester_conversion_depths = []

    def dfs_traverse(node, depth=0):
        nonlocal grignard_depths, ester_conversion_depths

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # Check for any Grignard reaction
            grignard_reaction_types = [
                "Grignard from ketone to alcohol",
                "Grignard from aldehyde to alcohol",
                "Grignard with CO2 to carboxylic acid",
                "Olefination of ketones with Grignard reagents",
                "Olefination of aldehydes with Grignard reagents",
                "Formation of Grignard reagents",
            ]

            # Check if any reactant contains "Mg" (indicating Grignard reagent)
            has_grignard_reagent = any("Mg" in r for r in reactants)

            if has_grignard_reagent or any(
                checker.check_reaction(rxn, rsmi) for rxn in grignard_reaction_types
            ):
                # Verify alcohol formation by checking product
                if (
                    checker.check_fg("Tertiary alcohol", product)
                    or (
                        checker.check_fg("Secondary alcohol", product)
                        and not any(
                            checker.check_fg("Secondary alcohol", r) for r in reactants
                        )
                    )
                    or (
                        checker.check_fg("Primary alcohol", product)
                        and not any(
                            checker.check_fg("Primary alcohol", r) for r in reactants
                        )
                    )
                ):
                    grignard_depths.append(depth)
                    print(
                        f"Grignard reaction forming alcohol found at depth {depth}, rsmi: {rsmi}"
                    )

            # Check for aryl halide or carboxylic acid to ester conversion
            if (
                (
                    any(checker.check_fg("Aromatic halide", r) for r in reactants)
                    or any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                )
                and checker.check_fg("Ester", product)
                and not any(checker.check_fg("Ester", r) for r in reactants)
            ):

                # Check for specific reactions that convert to esters
                ester_reaction_types = [
                    "Carbonylation with aryl formates",
                    "Oxidative esterification of primary alcohols",
                    "Esterification of Carboxylic Acids",
                    "Suzuki coupling",
                    "Heck reaction with vinyl ester",
                    "Oxidative Heck reaction with vinyl ester",
                    "Schotten-Baumann to ester",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Transesterification",
                    "Acetic anhydride and alcohol to ester",
                ]

                if any(
                    checker.check_reaction(rxn, rsmi) for rxn in ester_reaction_types
                ):
                    ester_conversion_depths.append(depth)
                    print(f"Conversion to ester found at depth {depth}, rsmi: {rsmi}")
                else:
                    # If no specific reaction type matches, but we have ester formation, still count it
                    ester_conversion_depths.append(depth)
                    print(
                        f"Generic ester formation found at depth {depth}, rsmi: {rsmi}"
                    )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if any Grignard happens early (higher depth) and any ester conversion happens late (lower depth)
    if grignard_depths and ester_conversion_depths:
        early_grignard = max(grignard_depths)  # Highest depth = earliest stage
        late_ester = min(ester_conversion_depths)  # Lowest depth = latest stage

        print(f"Found Grignard reactions at depths: {grignard_depths}")
        print(f"Found ester conversions at depths: {ester_conversion_depths}")

        if early_grignard > late_ester:
            print(
                f"Early Grignard (depth {early_grignard}), late ester conversion (depth {late_ester}) strategy detected"
            )
            return True
        else:
            print(
                f"Strategy not detected: Grignard at depth {early_grignard} is not earlier than ester conversion at depth {late_ester}"
            )
    else:
        if not grignard_depths:
            print("No Grignard reactions forming alcohols found")
        if not ester_conversion_depths:
            print("No conversions to esters found")

    return False
