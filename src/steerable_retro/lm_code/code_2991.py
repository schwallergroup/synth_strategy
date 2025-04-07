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
    Detects if the synthesis uses an early oxidation followed by later reductions.
    In retrosynthetic planning, early stage = high depth, late stage = low depth.
    """
    oxidation_depth = -1
    reduction_depths = []

    def dfs_traverse(node, current_depth=0):
        nonlocal oxidation_depth, reduction_depths

        # Try to get depth from metadata
        depth = current_depth
        if node.get("metadata", {}).get("depth", None) is not None:
            depth = node["metadata"]["depth"]
        elif "metadata" in node and "ID" in node["metadata"]:
            depth_info = str(node["metadata"]["ID"])
            if "Depth:" in depth_info:
                try:
                    depth = int(depth_info.split("Depth:")[1].split()[0])
                except Exception as e:
                    print(f"Error extracting depth from ID: {e}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Analyzing reaction at depth {depth}: {rsmi}")

            # In retrosynthesis, product is starting material, reactants are target products
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for oxidation reactions
                is_oxidation = False
                oxidation_reactions = [
                    "Oxidation of aldehydes to carboxylic acids",
                    "Oxidation of ketone to carboxylic acid",
                    "Oxidation of alcohol to carboxylic acid",
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                    "Oxidation of alkene to carboxylic acid",
                    "Oxidation of alkene to aldehyde",
                    "Oxidative esterification of primary alcohols",
                    "Oxidation of alcohol and aldehyde to ester",
                    "Quinone formation",
                    "Oxidation of boronic acids",
                    "Oxidation of boronic esters",
                ]

                for rxn_type in oxidation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_oxidation = True
                        print(f"Detected oxidation reaction '{rxn_type}' at depth {depth}")
                        break

                # Also check for alcohol oxidation by checking functional groups
                if not is_oxidation:
                    # In retrosynthesis, product is the starting material (reduced form)
                    # and reactants are the target products (oxidized form)
                    has_alcohol_in_product = (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                        or checker.check_fg("Aromatic alcohol", product)
                        or checker.check_fg("Phenol", product)
                    )

                    has_aldehyde_in_reactants = any(
                        checker.check_fg("Aldehyde", r) for r in reactants
                    )
                    has_acid_in_reactants = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    has_ketone_in_reactants = any(checker.check_fg("Ketone", r) for r in reactants)
                    has_ester_in_reactants = any(checker.check_fg("Ester", r) for r in reactants)

                    if has_alcohol_in_product and (
                        has_aldehyde_in_reactants
                        or has_acid_in_reactants
                        or has_ketone_in_reactants
                        or has_ester_in_reactants
                    ):
                        is_oxidation = True
                        print(f"Detected oxidation by functional group change at depth {depth}")

                if is_oxidation and (oxidation_depth == -1 or depth > oxidation_depth):
                    oxidation_depth = depth
                    print(f"Updated earliest oxidation depth to {depth}")

                # Check for reduction reactions
                is_reduction = False
                reduction_reactions = [
                    "Reduction of aldehydes and ketones to alcohols",
                    "Reduction of ester to primary alcohol",
                    "Reduction of ketone to secondary alcohol",
                    "Reduction of carboxylic acid to primary alcohol",
                    "Reduction of nitrile to amine",
                    "Reduction of primary amides to amines",
                    "Reduction of secondary amides to amines",
                    "Reduction of tertiary amides to amines",
                    "Reduction of nitro groups to amines",
                    "Hydrogenation (double to single)",
                    "Hydrogenation (triple to double)",
                    "Arene hydrogenation",
                    "Azide to amine reduction (Staudinger)",
                ]

                for rxn_type in reduction_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        is_reduction = True
                        print(f"Detected reduction reaction '{rxn_type}' at depth {depth}")
                        break

                # Also check for reduction by checking functional groups
                if not is_reduction:
                    # In retrosynthesis, product is the starting material (oxidized form)
                    # and reactants are the target products (reduced form)
                    has_carbonyl_in_product = (
                        checker.check_fg("Aldehyde", product)
                        or checker.check_fg("Ketone", product)
                        or checker.check_fg("Carboxylic acid", product)
                        or checker.check_fg("Ester", product)
                    )

                    has_alcohol_in_reactants = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        or checker.check_fg("Phenol", r)
                        for r in reactants
                    )

                    has_amine_in_reactants = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Tertiary amine", r)
                        for r in reactants
                    )

                    has_nitro_in_product = checker.check_fg("Nitro group", product)
                    has_nitrile_in_product = checker.check_fg("Nitrile", product)
                    has_azide_in_product = checker.check_fg("Azide", product)

                    has_alkene_in_product = checker.check_fg("Alkene", product) or checker.check_fg(
                        "Vinyl", product
                    )
                    has_alkyne_in_product = checker.check_fg("Alkyne", product)

                    # Check if reactants don't have the unsaturated groups (indicating reduction)
                    has_reduced_alkene = has_alkene_in_product and not any(
                        checker.check_fg("Alkene", r) or checker.check_fg("Vinyl", r)
                        for r in reactants
                    )

                    has_reduced_alkyne = has_alkyne_in_product and not any(
                        checker.check_fg("Alkyne", r) for r in reactants
                    )

                    if (
                        (has_carbonyl_in_product and has_alcohol_in_reactants)
                        or (has_nitro_in_product and has_amine_in_reactants)
                        or (
                            has_nitrile_in_product
                            and any(checker.check_fg("Primary amine", r) for r in reactants)
                        )
                        or (
                            has_azide_in_product
                            and any(checker.check_fg("Primary amine", r) for r in reactants)
                        )
                        or has_reduced_alkene
                        or has_reduced_alkyne
                    ):
                        is_reduction = True
                        print(f"Detected reduction by functional group change at depth {depth}")

                if is_reduction:
                    reduction_depths.append(depth)
                    print(f"Added reduction at depth {depth}")

            except Exception as e:
                print(f"Error processing SMILES in early_stage_redox_sequence: {e}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if oxidation occurs early and reduction occurs later
    if oxidation_depth != -1 and reduction_depths:
        # In retrosynthetic direction, higher depth = earlier in synthesis
        # So oxidation should have higher depth (earlier) than reduction (later)
        min_reduction_depth = min(reduction_depths)
        if oxidation_depth > min_reduction_depth:
            print(
                f"Detected early oxidation (depth {oxidation_depth}) followed by later reduction (depth {min_reduction_depth})"
            )
            return True
        else:
            print(
                f"Found oxidation and reduction, but oxidation (depth {oxidation_depth}) is not earlier than reduction (depth {min_reduction_depth})"
            )
    else:
        if oxidation_depth == -1:
            print("No oxidation reactions detected")
        if not reduction_depths:
            print("No reduction reactions detected")

    return False
