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
    Detects a synthesis route that includes both nitrile formation and reduction.
    """
    nitrile_formation = False
    nitrile_reduction = False

    # First check if any molecule in the route contains a nitrile group
    def check_for_nitrile(node):
        nonlocal nitrile_formation
        if node["type"] == "mol" and checker.check_fg("Nitrile", node["smiles"]):
            print(f"Nitrile found in molecule: {node['smiles']}")
            nitrile_formation = True
        for child in node.get("children", []):
            check_for_nitrile(child)

    check_for_nitrile(route)

    def dfs_traverse(node):
        nonlocal nitrile_formation, nitrile_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitrile formation
            if not nitrile_formation:
                # Check if product contains nitrile but reactants don't
                if checker.check_fg("Nitrile", product):
                    if not any(checker.check_fg("Nitrile", r) for r in reactants):
                        print(f"Nitrile formation detected (new nitrile in product): {rsmi}")
                        nitrile_formation = True

                # Check for specific nitrile formation reactions
                if checker.check_reaction("Schmidt reaction nitrile", rsmi):
                    print(f"Nitrile formation reaction detected (Schmidt reaction): {rsmi}")
                    nitrile_formation = True

                # Check for reverse reactions where nitrile is consumed
                # These indicate nitrile was formed earlier in the synthesis (forward direction)
                if any(checker.check_fg("Nitrile", r) for r in reactants):
                    if checker.check_reaction(
                        "Oxidation of nitrile to carboxylic acid", rsmi
                    ) or checker.check_reaction("Nitrile to amide", rsmi):
                        print(f"Nitrile formation inferred from consumption reaction: {rsmi}")
                        nitrile_formation = True

            # Check for nitrile reduction
            if not nitrile_reduction:
                # Check if reactants contain nitrile and product contains amine
                if any(checker.check_fg("Nitrile", r) for r in reactants) and (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                ):
                    print(f"Potential nitrile reduction detected in reaction: {rsmi}")
                    nitrile_reduction = True

                # Check for specific nitrile reduction reactions
                if checker.check_reaction("Reduction of nitrile to amine", rsmi):
                    print(f"Nitrile reduction reaction detected: {rsmi}")
                    nitrile_reduction = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Nitrile formation: {nitrile_formation}")
    print(f"Nitrile reduction: {nitrile_reduction}")

    # Check if both conditions are met
    return nitrile_formation and nitrile_reduction
