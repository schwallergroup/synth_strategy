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
    Detects a linear synthesis strategy with early ring formation, mid-stage cross-coupling,
    and late-stage N-alkylation.
    """
    # Initialize flags for each key transformation
    found_ring_formation = False
    found_nitrile_to_acid = False
    found_acid_to_ester = False
    found_borylation = False
    found_suzuki_coupling = False
    found_n_alkylation = False

    # Track depths for ordering
    ring_formation_depth = -1
    nitrile_to_acid_depth = -1
    acid_to_ester_depth = -1
    borylation_depth = -1
    suzuki_depth = -1
    n_alkylation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_ring_formation, found_nitrile_to_acid, found_acid_to_ester
        nonlocal found_borylation, found_suzuki_coupling, found_n_alkylation
        nonlocal ring_formation_depth, nitrile_to_acid_depth, acid_to_ester_depth
        nonlocal borylation_depth, suzuki_depth, n_alkylation_depth

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}: {rsmi}")

                # Check for ring formation (early stage)
                if depth >= 3:  # High depth = early in synthesis (relaxed constraint)
                    # Check for common rings that might be formed
                    rings_to_check = [
                        "benzene",
                        "pyridine",
                        "furan",
                        "thiophene",
                        "pyrrole",
                        "imidazole",
                        "oxazole",
                        "thiazole",
                        "pyrimidine",
                        "cyclopentane",
                        "cyclohexane",
                    ]

                    # Check if product has a ring
                    for ring in rings_to_check:
                        if checker.check_ring(ring, product_smiles):
                            found_ring_formation = True
                            ring_formation_depth = depth
                            print(f"Detected {ring} in product at depth {depth}")
                            break

                # Check for nitrile to acid conversion
                if any(
                    checker.check_fg("Nitrile", r) for r in reactants_smiles
                ) and checker.check_fg("Carboxylic acid", product_smiles):
                    found_nitrile_to_acid = True
                    nitrile_to_acid_depth = depth
                    print(
                        f"Detected nitrile to carboxylic acid conversion at depth {depth}"
                    )

                # Check for acid to ester conversion
                if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                    found_acid_to_ester = True
                    acid_to_ester_depth = depth
                    print(
                        f"Detected esterification of carboxylic acids at depth {depth}"
                    )
                elif any(
                    checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                ) and checker.check_fg("Ester", product_smiles):
                    found_acid_to_ester = True
                    acid_to_ester_depth = depth
                    print(f"Detected acid to ester conversion at depth {depth}")

                # Check for borylation (expanded to include boronic esters)
                if any(
                    checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                ) and (
                    checker.check_fg("Boronic acid", product_smiles)
                    or checker.check_fg("Boronic ester", product_smiles)
                ):
                    found_borylation = True
                    borylation_depth = depth
                    print(f"Detected borylation at depth {depth}")

                # Check for Suzuki coupling
                suzuki_reactions = [
                    "Suzuki coupling with boronic acids",
                    "Suzuki coupling with boronic acids OTf",
                    "Suzuki coupling with sulfonic esters",
                    "Suzuki coupling with boronic esters",
                    "Suzuki coupling with boronic esters OTf",
                    "Suzuki",  # Added generic Suzuki reaction
                ]

                for rxn_type in suzuki_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        found_suzuki_coupling = True
                        suzuki_depth = depth
                        print(f"Detected {rxn_type} at depth {depth}")
                        break

                # If no specific Suzuki reaction was detected, check for patterns
                if not found_suzuki_coupling and depth < 6:  # Mid-stage
                    # Check if reactants contain boronic acid/ester and aryl halide
                    has_boronic = any(
                        checker.check_fg("Boronic acid", r)
                        or checker.check_fg("Boronic ester", r)
                        for r in reactants_smiles
                    )
                    has_aryl_halide = any(
                        checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                    )

                    if has_boronic and has_aryl_halide:
                        found_suzuki_coupling = True
                        suzuki_depth = depth
                        print(
                            f"Detected likely Suzuki coupling pattern at depth {depth}"
                        )

                # Check for N-alkylation (late stage) - relaxed depth constraint
                if depth <= 3:  # Low depth = late in synthesis
                    n_alkylation_reactions = [
                        "N-alkylation of primary amines with alkyl halides",
                        "N-alkylation of secondary amines with alkyl halides",
                        "Alkylation of amines",
                        "N-alkylation",  # Added generic N-alkylation
                    ]

                    for rxn_type in n_alkylation_reactions:
                        if checker.check_reaction(rxn_type, rsmi):
                            found_n_alkylation = True
                            n_alkylation_depth = depth
                            print(f"Detected {rxn_type} at depth {depth}")
                            break

                    # If no specific N-alkylation reaction was detected, check for patterns
                    if not found_n_alkylation:
                        # Check if product has N-alkyl group that wasn't in reactants
                        if checker.check_fg(
                            "Secondary amine", product_smiles
                        ) or checker.check_fg("Tertiary amine", product_smiles):
                            if any(
                                checker.check_fg("Primary amine", r)
                                or checker.check_fg("Secondary amine", r)
                                for r in reactants_smiles
                            ):
                                found_n_alkylation = True
                                n_alkylation_depth = depth
                                print(f"Detected N-alkylation pattern at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all transformations were found with a more flexible ordering
    # We still want ring formation early, N-alkylation late, and the rest in between
    correct_strategy = (
        found_ring_formation
        and found_nitrile_to_acid
        and found_acid_to_ester
        and (found_borylation or found_suzuki_coupling)
        and found_n_alkylation  # Allow either borylation or Suzuki
    )

    # Check relative ordering with some flexibility
    correct_ordering = (
        (ring_formation_depth > n_alkylation_depth)
        and (  # Ring formation before N-alkylation
            nitrile_to_acid_depth > n_alkylation_depth
        )
        and (  # Nitrile to acid before N-alkylation
            acid_to_ester_depth > n_alkylation_depth
        )  # Acid to ester before N-alkylation
    )

    print(f"Strategy detected: {correct_strategy and correct_ordering}")
    print(f"Ring formation: {found_ring_formation} at depth {ring_formation_depth}")
    print(f"Nitrile to acid: {found_nitrile_to_acid} at depth {nitrile_to_acid_depth}")
    print(f"Acid to ester: {found_acid_to_ester} at depth {acid_to_ester_depth}")
    print(f"Borylation: {found_borylation} at depth {borylation_depth}")
    print(f"Suzuki coupling: {found_suzuki_coupling} at depth {suzuki_depth}")
    print(f"N-alkylation: {found_n_alkylation} at depth {n_alkylation_depth}")

    return correct_strategy and correct_ordering
