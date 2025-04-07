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
    This function detects if the synthesis involves a specific sequence of
    heteroatom bond formations: C-O → C-B → C-C → C-N → C-N
    """
    # Define patterns for each bond type
    co_pattern = "C-O formation"  # C-O bond formation
    cb_pattern = "C-B formation"  # Borylation
    cc_pattern = "C-C formation"  # C-C bond formation
    cn1_pattern = "C-N formation 1"  # First C-N bond formation
    cn2_pattern = "C-N formation 2"  # Second C-N bond formation

    # Track bond formations in order (from early to late)
    bond_formation_sequence = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Extract reactants and product for manual analysis
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for C-O bond formation
                co_formation = False
                if (
                    checker.check_reaction("Esterification of Carboxylic Acids", rsmi)
                    or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                    or checker.check_reaction("Alcohol protection with silyl ethers", rsmi)
                    or checker.check_reaction(
                        "O-alkylation of carboxylic acids with diazo compounds", rsmi
                    )
                    or checker.check_reaction("Oxidative esterification of primary alcohols", rsmi)
                    or checker.check_reaction("Acetic anhydride and alcohol to ester", rsmi)
                    or checker.check_reaction("Mitsunobu esterification", rsmi)
                    or checker.check_reaction("Mitsunobu aryl ether", rsmi)
                ):
                    co_formation = True

                # Manual check for C-O bond formation (esterification)
                if not co_formation:
                    # Check for carboxylic acid and alcohol to ester
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    has_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                    )
                    has_ester = checker.check_fg("Ester", product)

                    if has_carboxylic_acid and has_alcohol and has_ester:
                        co_formation = True

                if co_formation:
                    bond_formation_sequence.append((depth, co_pattern))
                    print(f"C-O bond formation detected at depth {depth}")

                # Check for C-B bond formation
                cb_formation = False
                if (
                    checker.check_reaction("Preparation of boronic acids", rsmi)
                    or checker.check_reaction(
                        "Preparation of boronic acids without boronic ether", rsmi
                    )
                    or checker.check_reaction(
                        "Preparation of boronic acids from trifluoroborates", rsmi
                    )
                    or checker.check_reaction("Preparation of boronic esters", rsmi)
                    or checker.check_reaction("Synthesis of boronic acids", rsmi)
                ):
                    cb_formation = True

                # Manual check for C-B bond formation
                if not cb_formation:
                    has_boronic_acid_product = checker.check_fg("Boronic acid", product)
                    has_boronic_ester_product = checker.check_fg("Boronic ester", product)

                    if has_boronic_acid_product or has_boronic_ester_product:
                        cb_formation = True

                if cb_formation:
                    bond_formation_sequence.append((depth, cb_pattern))
                    print(f"C-B bond formation detected at depth {depth}")

                # Check for C-C bond formation
                cc_formation = False
                if (
                    checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                    or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
                    or checker.check_reaction("Negishi coupling", rsmi)
                    or checker.check_reaction("Heck terminal vinyl", rsmi)
                    or checker.check_reaction("Stille reaction_aryl", rsmi)
                    or checker.check_reaction("Grignard from aldehyde to alcohol", rsmi)
                    or checker.check_reaction("Grignard from ketone to alcohol", rsmi)
                    or checker.check_reaction("Wittig reaction with triphenylphosphorane", rsmi)
                    or checker.check_reaction("Aldol condensation", rsmi)
                    or checker.check_reaction("Suzuki", rsmi)
                ):
                    cc_formation = True

                # Manual check for C-C bond formation (Suzuki coupling)
                if not cc_formation:
                    has_boronic_reactant = any(
                        checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                        for r in reactants
                    )
                    has_aryl_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

                    if has_boronic_reactant and has_aryl_halide:
                        cc_formation = True

                if cc_formation:
                    bond_formation_sequence.append((depth, cc_pattern))
                    print(f"C-C bond formation detected at depth {depth}")

                # Check for first C-N bond formation
                cn1_formation = False
                if (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                    or checker.check_reaction(
                        "N-alkylation of primary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "N-alkylation of secondary amines with alkyl halides", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation primary amine", rsmi
                    )
                    or checker.check_reaction(
                        "Buchwald-Hartwig/Ullmann-Goldberg/N-arylation secondary amine", rsmi
                    )
                    or checker.check_reaction(
                        "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)", rsmi
                    )
                    or checker.check_reaction("Alkylation of amines", rsmi)
                    or checker.check_reaction("reductive amination", rsmi)
                ):
                    cn1_formation = True

                # Manual check for C-N bond formation (reductive amination)
                if not cn1_formation:
                    has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
                    has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                    has_secondary_amine_product = checker.check_fg("Secondary amine", product)

                    if has_aldehyde and has_primary_amine and has_secondary_amine_product:
                        cn1_formation = True

                if cn1_formation:
                    bond_formation_sequence.append((depth, cn1_pattern))
                    print(f"First C-N bond formation detected at depth {depth}")

                # Check for second C-N bond formation (amide formation)
                cn2_formation = False
                if (
                    checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                    or checker.check_reaction("Acylation of primary amines", rsmi)
                    or checker.check_reaction("Acylation of secondary amines", rsmi)
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                ):
                    cn2_formation = True

                # Manual check for C-N bond formation (amide formation)
                if not cn2_formation:
                    has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )
                    has_amide = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    if (has_acyl_halide or has_amine) and has_amide:
                        cn2_formation = True

                if cn2_formation:
                    bond_formation_sequence.append((depth, cn2_pattern))
                    print(f"Second C-N bond formation (amide) detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Sort by depth (high to low, as high depth is early in synthesis)
    bond_formation_sequence.sort(key=lambda x: x[0], reverse=True)
    print(f"Bond formation sequence (sorted by depth): {bond_formation_sequence}")

    # Extract just the patterns in sequence
    patterns = [p for _, p in bond_formation_sequence]
    print(f"Patterns in sequence: {patterns}")

    # Check if the sequence matches our expected pattern
    expected_sequence = [co_pattern, cb_pattern, cc_pattern, cn1_pattern, cn2_pattern]
    print(f"Expected sequence: {expected_sequence}")

    # Check if the patterns appear in the correct order (not necessarily consecutive)
    last_found_idx = -1
    sequence_valid = True

    # Track which patterns we've found
    found_expected_patterns = []

    for pattern in expected_sequence:
        if pattern in patterns:
            current_idx = patterns.index(pattern)
            found_expected_patterns.append(pattern)
            if current_idx <= last_found_idx:
                print(
                    f"Heteroatom bond formation sequence not in expected order: {pattern} found at position {current_idx}, after position {last_found_idx}"
                )
                sequence_valid = False
                break
            last_found_idx = current_idx

    # Check if we found at least 3 of the expected patterns in the correct order
    found_patterns = len(found_expected_patterns)
    print(f"Found {found_patterns} patterns from the expected sequence: {found_expected_patterns}")

    if found_patterns >= 3 and sequence_valid:
        print(
            f"Sequential heteroatom bond formation detected with {found_patterns} steps in the expected order"
        )
        return True

    print(f"Insufficient heteroatom bond formations in the expected sequence")
    return False
