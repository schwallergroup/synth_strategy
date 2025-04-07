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
    Detects a synthesis strategy involving a methoxybenzyl-indazole core
    with a specific functional group transformation sequence:
    aldehyde → vinyl → diol → mesylate
    """
    # Track the functional group transformations with their depths
    transformations = []
    has_methoxybenzyl_indazole = False

    def dfs_traverse(node, depth=0):
        nonlocal transformations, has_methoxybenzyl_indazole

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check for methoxybenzyl-indazole core in any molecule
            if (
                checker.check_ring("indazole", mol_smiles)
                and "OC" in mol_smiles
                and checker.check_fg("Ether", mol_smiles)
            ):
                has_methoxybenzyl_indazole = True
                print(f"Detected methoxybenzyl-indazole core in molecule: {mol_smiles}")

        elif (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for aldehyde → vinyl transformation
            if (
                any(checker.check_fg("Aldehyde", r) for r in reactants)
                and checker.check_fg("Vinyl", product)
                and not checker.check_fg("Aldehyde", product)
            ):
                transformations.append(("aldehyde_to_vinyl", depth))
                print(
                    f"Detected aldehyde to vinyl transformation at depth {depth}: {rsmi}"
                )

            # Check for vinyl → diol transformation
            if (
                any(checker.check_fg("Vinyl", r) for r in reactants)
                and not checker.check_fg("Vinyl", product)
                and (
                    (
                        checker.check_fg("Primary alcohol", product)
                        and checker.check_fg("Secondary alcohol", product)
                    )
                    or product.count("OH") >= 2  # Backup check for diols
                )
            ):
                transformations.append(("vinyl_to_diol", depth))
                print(f"Detected vinyl to diol transformation at depth {depth}: {rsmi}")

            # Check for diol → mesylate transformation
            if any(
                r.count("OH") >= 2
                or (
                    checker.check_fg("Primary alcohol", r)
                    and checker.check_fg("Secondary alcohol", r)
                )
                for r in reactants
            ) and checker.check_fg("Mesylate", product):
                transformations.append(("diol_to_mesylate", depth))
                print(
                    f"Detected diol to mesylate transformation at depth {depth}: {rsmi}"
                )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if all required transformations are present
    required_transformations = [
        "aldehyde_to_vinyl",
        "vinyl_to_diol",
        "diol_to_mesylate",
    ]
    detected_transformations = [t[0] for t in transformations]
    all_present = all(t in detected_transformations for t in required_transformations)

    # In retrosynthesis, the order should be reversed (higher depth = earlier stage)
    # Create a dictionary mapping transformation types to their depths
    transformation_depths = {t[0]: t[1] for t in transformations}

    # Check if the transformations appear in the correct order for retrosynthesis
    # In retrosynthesis: mesylate (late stage, low depth) -> diol -> vinyl -> aldehyde (early stage, high depth)
    correct_order = True
    if all_present:
        # Check that diol_to_mesylate occurs at a lower depth (later stage) than vinyl_to_diol
        if (
            transformation_depths["diol_to_mesylate"]
            >= transformation_depths["vinyl_to_diol"]
        ):
            correct_order = False
            print(
                f"Incorrect order: diol_to_mesylate (depth {transformation_depths['diol_to_mesylate']}) should occur at lower depth than vinyl_to_diol (depth {transformation_depths['vinyl_to_diol']})"
            )

        # Check that vinyl_to_diol occurs at a lower depth (later stage) than aldehyde_to_vinyl
        if (
            transformation_depths["vinyl_to_diol"]
            >= transformation_depths["aldehyde_to_vinyl"]
        ):
            correct_order = False
            print(
                f"Incorrect order: vinyl_to_diol (depth {transformation_depths['vinyl_to_diol']}) should occur at lower depth than aldehyde_to_vinyl (depth {transformation_depths['aldehyde_to_vinyl']})"
            )
    else:
        correct_order = False
        missing = [
            t for t in required_transformations if t not in detected_transformations
        ]
        print(f"Missing transformations: {missing}")

    strategy_detected = has_methoxybenzyl_indazole and all_present and correct_order

    if strategy_detected:
        print("Detected complete functional group sequence strategy")
    else:
        print(
            f"Strategy detection failed: has_core={has_methoxybenzyl_indazole}, all_present={all_present}, correct_order={correct_order}"
        )
        print(f"Detected transformations with depths: {transformations}")

    return strategy_detected
