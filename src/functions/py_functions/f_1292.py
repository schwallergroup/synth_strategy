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
    This function detects a linear synthesis strategy with a sequence of
    functional group interconversions, particularly involving carboxylic acids,
    esters, and other carbonyl derivatives.
    """
    # Track functional group transformations
    transformations = []
    linear_synthesis = True

    def dfs_traverse(node, depth=0):
        nonlocal transformations, linear_synthesis

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants_smiles = parts[0].split(".")
            product_smiles = parts[2]

            # Check if synthesis is linear (typically 1-2 reactants per step)
            if len(reactants_smiles) > 2:
                linear_synthesis = False

            # Identify functional groups in reactants and products
            reactant_groups = []
            product_groups = []
            reaction_type = None

            # Check reactants
            for reactant in reactants_smiles:
                try:
                    if checker.check_fg("Carboxylic acid", reactant):
                        reactant_groups.append("carboxylic_acid")
                    if checker.check_fg("Ester", reactant):
                        reactant_groups.append("ester")
                    if checker.check_fg("Acyl halide", reactant):
                        reactant_groups.append("acid_chloride")
                    if checker.check_fg("Enol", reactant):
                        reactant_groups.append("enol")
                    if checker.check_fg("Anhydride", reactant):
                        reactant_groups.append("anhydride")
                    if (
                        checker.check_fg("Primary amide", reactant)
                        or checker.check_fg("Secondary amide", reactant)
                        or checker.check_fg("Tertiary amide", reactant)
                    ):
                        reactant_groups.append("amide")
                except Exception as e:
                    print(f"Error checking reactant functional groups: {e}")
                    continue

            # Check product
            try:
                if checker.check_fg("Carboxylic acid", product_smiles):
                    product_groups.append("carboxylic_acid")
                if checker.check_fg("Ester", product_smiles):
                    product_groups.append("ester")
                if checker.check_fg("Acyl halide", product_smiles):
                    product_groups.append("acid_chloride")
                if checker.check_fg("Enol", product_smiles):
                    product_groups.append("enol")
                if checker.check_fg("Anhydride", product_smiles):
                    product_groups.append("anhydride")
                if (
                    checker.check_fg("Primary amide", product_smiles)
                    or checker.check_fg("Secondary amide", product_smiles)
                    or checker.check_fg("Tertiary amide", product_smiles)
                ):
                    product_groups.append("amide")
            except Exception as e:
                print(f"Error checking product functional groups: {e}")
                pass

            # Identify reaction type
            if checker.check_reaction("Esterification of Carboxylic Acids", rsmi):
                reaction_type = "esterification"
            elif checker.check_reaction("Schotten-Baumann to ester", rsmi):
                reaction_type = "schotten_baumann_ester"
            elif checker.check_reaction(
                "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi
            ):
                reaction_type = "ester_hydrolysis"
            elif checker.check_reaction(
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                rsmi,
            ):
                reaction_type = "acylation"
            elif checker.check_reaction(
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                rsmi,
            ):
                reaction_type = "acylation"
            elif checker.check_reaction("Transesterification", rsmi):
                reaction_type = "transesterification"
            elif checker.check_reaction(
                "Hydrogenolysis of amides/imides/carbamates", rsmi
            ) or checker.check_reaction("Hydrolysis of amides/imides/carbamates", rsmi):
                reaction_type = "amide_hydrolysis"
            elif (
                checker.check_reaction("Acyl chlorides from alcohols", rsmi)
                or checker.check_reaction("Acyl bromides from alcohols", rsmi)
                or checker.check_reaction("Acyl iodides from alcohols", rsmi)
            ):
                reaction_type = "acyl_halide_formation"
            elif checker.check_reaction(
                "Carboxylic acid with primary amine to amide", rsmi
            ):
                reaction_type = "amide_formation"
            elif (
                checker.check_reaction("Ester with primary amine to amide", rsmi)
                or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                or checker.check_reaction("Ester with ammonia to amide", rsmi)
            ):
                reaction_type = "ester_to_amide"

            # Record transformation if we have functional groups in both reactants and products
            if reactant_groups and product_groups:
                transformation = {
                    "depth": depth,
                    "from": reactant_groups,
                    "to": product_groups,
                    "reaction_type": reaction_type,
                    "rsmi": rsmi,
                }
                transformations.append(transformation)
                print(
                    f"Detected transformation at depth {depth}: {reactant_groups} → {product_groups} (Reaction: {reaction_type})"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check for functional group interconversion sequence
    has_functional_group_sequence = False
    if len(transformations) >= 2 and linear_synthesis:
        # Sort transformations by depth (higher depth = earlier in synthesis)
        sorted_transformations = sorted(
            transformations, key=lambda x: x["depth"], reverse=True
        )

        # Check for patterns of functional group interconversions
        pattern_count = 0
        for i in range(len(sorted_transformations)):
            current = sorted_transformations[i]

            # Check for acid → ester
            if "carboxylic_acid" in current["from"] and "ester" in current["to"]:
                pattern_count += 1
                print(f"Found acid → ester pattern at depth {current['depth']}")

            # Check for acid chloride → ester
            elif "acid_chloride" in current["from"] and "ester" in current["to"]:
                pattern_count += 1
                print(
                    f"Found acid chloride → ester pattern at depth {current['depth']}"
                )

            # Check for ester → acid
            elif "ester" in current["from"] and "carboxylic_acid" in current["to"]:
                pattern_count += 1
                print(f"Found ester → acid pattern at depth {current['depth']}")

            # Check for ester → ester (transesterification)
            elif (
                "ester" in current["from"]
                and "ester" in current["to"]
                and len(current["from"]) == 1
            ):
                pattern_count += 1
                print(f"Found ester → ester pattern at depth {current['depth']}")

            # Check for acid → acid chloride
            elif (
                "carboxylic_acid" in current["from"]
                and "acid_chloride" in current["to"]
            ):
                pattern_count += 1
                print(f"Found acid → acid chloride pattern at depth {current['depth']}")

            # Check for acid → amide
            elif "carboxylic_acid" in current["from"] and "amide" in current["to"]:
                pattern_count += 1
                print(f"Found acid → amide pattern at depth {current['depth']}")

            # Check for acid chloride → amide
            elif "acid_chloride" in current["from"] and "amide" in current["to"]:
                pattern_count += 1
                print(
                    f"Found acid chloride → amide pattern at depth {current['depth']}"
                )

            # Check for amide → acid
            elif "amide" in current["from"] and "carboxylic_acid" in current["to"]:
                pattern_count += 1
                print(f"Found amide → acid pattern at depth {current['depth']}")

            # Check for ester → amide
            elif "ester" in current["from"] and "amide" in current["to"]:
                pattern_count += 1
                print(f"Found ester → amide pattern at depth {current['depth']}")

            # Check for anhydride → acid
            elif "anhydride" in current["from"] and "carboxylic_acid" in current["to"]:
                pattern_count += 1
                print(f"Found anhydride → acid pattern at depth {current['depth']}")

            # Check for anhydride → ester
            elif "anhydride" in current["from"] and "ester" in current["to"]:
                pattern_count += 1
                print(f"Found anhydride → ester pattern at depth {current['depth']}")

            # Check for anhydride → amide
            elif "anhydride" in current["from"] and "amide" in current["to"]:
                pattern_count += 1
                print(f"Found anhydride → amide pattern at depth {current['depth']}")

        # We need at least 2 functional group interconversions in sequence
        has_functional_group_sequence = pattern_count >= 2

        # Check if the interconversions form a sequence
        # This means at least one product functional group should be used as a reactant in a later step
        if pattern_count >= 2:
            has_sequence = False
            for i in range(len(sorted_transformations) - 1):
                current = sorted_transformations[i]
                for j in range(i + 1, len(sorted_transformations)):
                    next_step = sorted_transformations[j]
                    # Check if any product functional group from current step is used as reactant in next step
                    for fg in current["to"]:
                        if fg in next_step["from"]:
                            has_sequence = True
                            print(
                                f"Found sequence: {fg} from depth {current['depth']} used at depth {next_step['depth']}"
                            )
                            break
                    if has_sequence:
                        break
                if has_sequence:
                    break

            has_functional_group_sequence = has_sequence

    strategy_detected = linear_synthesis and has_functional_group_sequence

    if strategy_detected:
        print(
            f"Detected linear synthesis with functional group interconversion sequence"
        )
        for t in sorted(transformations, key=lambda x: x["depth"], reverse=True):
            print(
                f"  Depth {t['depth']}: {t['from']} → {t['to']} (Reaction: {t['reaction_type']})"
            )
    else:
        print(f"Linear synthesis with functional group interconversion NOT detected")
        print(f"  Linear synthesis: {linear_synthesis}")
        print(f"  Has functional group sequence: {has_functional_group_sequence}")
        print(f"  Number of transformations: {len(transformations)}")

    return strategy_detected
