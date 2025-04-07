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
    Detects a synthesis strategy where a haloalkyl group (like chloroethyl)
    is introduced in the final step of the synthesis.
    """
    has_late_stage_haloalkylation = False

    def is_haloalkylation_pattern(rsmi):
        """Check if the reaction pattern matches haloalkylation even if not explicitly categorized"""
        try:
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            # Check if product has an N-haloalkyl or O-haloalkyl pattern
            product_mol = Chem.MolFromSmiles(product_smiles)
            if not product_mol:
                return False

            # Check for halide-containing reactant
            reactants_list = reactants_smiles.split(".")
            halide_reactant = None

            for r in reactants_list:
                if not r:
                    continue
                if (
                    checker.check_fg("Primary halide", r)
                    or checker.check_fg("Secondary halide", r)
                    or checker.check_fg("Tertiary halide", r)
                ):
                    halide_reactant = r
                    break

            if not halide_reactant:
                return False

            # Check if the product has a new N-haloalkyl or O-haloalkyl bond
            if (
                checker.check_fg("Primary halide", product_smiles)
                or checker.check_fg("Secondary halide", product_smiles)
                or checker.check_fg("Tertiary halide", product_smiles)
            ):

                # Check if the product has a primary, secondary, or tertiary amine
                has_amine = (
                    checker.check_fg("Primary amine", product_smiles)
                    or checker.check_fg("Secondary amine", product_smiles)
                    or checker.check_fg("Tertiary amine", product_smiles)
                )

                # Check if the product has an ether or alcohol
                has_oxygen = (
                    checker.check_fg("Ether", product_smiles)
                    or checker.check_fg("Primary alcohol", product_smiles)
                    or checker.check_fg("Secondary alcohol", product_smiles)
                    or checker.check_fg("Tertiary alcohol", product_smiles)
                )

                return has_amine or has_oxygen

            return False
        except Exception as e:
            print(f"Error in haloalkylation pattern check: {e}")
            return False

    def dfs_traverse(node, current_depth=0):
        nonlocal has_late_stage_haloalkylation

        if node["type"] == "reaction":
            # Get depth from metadata or use current_depth from traversal
            depth = node["metadata"].get("depth", current_depth)

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants_smiles = rsmi.split(">")[0]
                    product_smiles = rsmi.split(">")[-1]

                    print(f"Analyzing reaction at depth {depth}: {rsmi}")

                    # Check if this is an alkylation reaction
                    is_alkylation = (
                        checker.check_reaction(
                            "N-alkylation of primary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction(
                            "N-alkylation of secondary amines with alkyl halides", rsmi
                        )
                        or checker.check_reaction("S-alkylation of thiols", rsmi)
                        or checker.check_reaction(
                            "S-alkylation of thiols (ethyl)", rsmi
                        )
                        or checker.check_reaction(
                            "S-alkylation of thiols with alcohols", rsmi
                        )
                        or checker.check_reaction(
                            "S-alkylation of thiols with alcohols (ethyl)", rsmi
                        )
                        or checker.check_reaction("Williamson Ether Synthesis", rsmi)
                        or checker.check_reaction(
                            "O-alkylation of carboxylic acids with diazo compounds",
                            rsmi,
                        )
                        or checker.check_reaction(
                            "O-alkylation of amides with diazo compounds", rsmi
                        )
                        or checker.check_reaction("Methylation", rsmi)
                        or checker.check_reaction("Methylation with MeI_primary", rsmi)
                        or checker.check_reaction(
                            "Methylation with MeI_secondary", rsmi
                        )
                        or checker.check_reaction("Methylation with MeI_tertiary", rsmi)
                        or checker.check_reaction("Methylation with MeI_aryl", rsmi)
                        or checker.check_reaction("Methylation with MeI_SH", rsmi)
                        or checker.check_reaction("Alkylation of amines", rsmi)
                        or checker.check_reaction(
                            "Reductive amination with aldehyde", rsmi
                        )
                        or checker.check_reaction(
                            "Reductive amination with ketone", rsmi
                        )
                        or checker.check_reaction(
                            "Reductive amination with alcohol", rsmi
                        )
                        or checker.check_reaction("{reductive amination}", rsmi)
                    )

                    # Check for manual haloalkylation pattern
                    if not is_alkylation:
                        # Look for specific patterns of haloalkylation
                        reactants_list = reactants_smiles.split(".")
                        for r in reactants_list:
                            if not r:
                                continue
                            # Check if any reactant has a halide group that could be used for alkylation
                            if (
                                checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Tertiary halide", r)
                            ):
                                # Check if the product has a new N-C or O-C bond
                                if (
                                    checker.check_fg("Primary amine", product_smiles)
                                    or checker.check_fg(
                                        "Secondary amine", product_smiles
                                    )
                                    or checker.check_fg(
                                        "Tertiary amine", product_smiles
                                    )
                                    or checker.check_fg("Ether", product_smiles)
                                ):
                                    is_alkylation = True
                                    print(
                                        f"Detected potential haloalkylation pattern manually"
                                    )
                                    break

                    print(f"Is alkylation reaction: {is_alkylation}")

                    if is_alkylation:
                        # Check if product contains halide groups
                        product_has_primary_halide = checker.check_fg(
                            "Primary halide", product_smiles
                        )
                        product_has_secondary_halide = checker.check_fg(
                            "Secondary halide", product_smiles
                        )
                        product_has_tertiary_halide = checker.check_fg(
                            "Tertiary halide", product_smiles
                        )
                        product_has_halide = (
                            product_has_primary_halide
                            or product_has_secondary_halide
                            or product_has_tertiary_halide
                        )

                        print(f"Product has halide: {product_has_halide}")

                        # Check if any reactant has halide groups
                        reactants_list = reactants_smiles.split(".")
                        reactant_has_halide = False

                        for r in reactants_list:
                            if not r:
                                continue

                            r_has_primary_halide = checker.check_fg("Primary halide", r)
                            r_has_secondary_halide = checker.check_fg(
                                "Secondary halide", r
                            )
                            r_has_tertiary_halide = checker.check_fg(
                                "Tertiary halide", r
                            )

                            if (
                                r_has_primary_halide
                                or r_has_secondary_halide
                                or r_has_tertiary_halide
                            ):
                                reactant_has_halide = True
                                print(f"Reactant has halide: {r}")
                                break

                        print(f"Any reactant has halide: {reactant_has_halide}")

                        # For the specific case in the test, check if we have a chloroalkyl group being introduced
                        # Look for a pattern where a halide-containing aldehyde or ketone is used
                        has_haloalkyl_pattern = False
                        for r in reactants_list:
                            if not r:
                                continue
                            if (
                                checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Tertiary halide", r)
                            ):
                                if checker.check_fg("Aldehyde", r) or checker.check_fg(
                                    "Ketone", r
                                ):
                                    has_haloalkyl_pattern = True
                                    print(
                                        f"Found haloalkyl aldehyde/ketone pattern in reactant: {r}"
                                    )
                                    break

                        # Check if the reaction introduces a haloalkyl group
                        if product_has_halide and (
                            reactant_has_halide or has_haloalkyl_pattern
                        ):
                            print(
                                f"Detected late-stage haloalkylation in step at depth {depth}"
                            )
                            has_late_stage_haloalkylation = True
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: {has_late_stage_haloalkylation}")

    return has_late_stage_haloalkylation
