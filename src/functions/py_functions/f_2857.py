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
    This function detects if the synthesis route involves nitrile chemistry,
    particularly the introduction of a nitrile group or transformations of nitriles.
    """
    nitrile_chemistry_detected = False

    def dfs_traverse(node):
        nonlocal nitrile_chemistry_detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for nitrile introduction
            product_has_nitrile = checker.check_fg("Nitrile", product)

            if product_has_nitrile:
                # Check if nitrile is introduced (not present in all reactants)
                all_reactants_have_nitrile = all(
                    checker.check_fg("Nitrile", r) for r in reactants
                )

                if not all_reactants_have_nitrile:
                    print(f"Nitrile introduction detected in reaction: {rsmi}")
                    nitrile_chemistry_detected = True

            # Check for nitrile transformations
            reactants_have_nitrile = any(
                checker.check_fg("Nitrile", r) for r in reactants
            )

            if reactants_have_nitrile:
                # Check for common nitrile transformations
                if checker.check_reaction("Nitrile to amide", rsmi):
                    print(f"Nitrile to amide transformation detected: {rsmi}")
                    nitrile_chemistry_detected = True
                elif checker.check_reaction(
                    "Oxidation of nitrile to carboxylic acid", rsmi
                ):
                    print(f"Nitrile to carboxylic acid transformation detected: {rsmi}")
                    nitrile_chemistry_detected = True
                elif checker.check_reaction("Reduction of nitrile to amine", rsmi):
                    print(f"Nitrile to amine reduction detected: {rsmi}")
                    nitrile_chemistry_detected = True
                elif checker.check_reaction("Grignard from nitrile to ketone", rsmi):
                    print(f"Nitrile to ketone via Grignard detected: {rsmi}")
                    nitrile_chemistry_detected = True

                # Check if nitrile is being used in heterocycle formation
                if not product_has_nitrile and (
                    checker.check_reaction(
                        "Azide-nitrile click cycloaddition to tetrazole", rsmi
                    )
                    or checker.check_reaction(
                        "Azide-nitrile click cycloaddition to triazole", rsmi
                    )
                    or checker.check_reaction("tetrazole_terminal", rsmi)
                ):
                    print(f"Nitrile used in heterocycle formation: {rsmi}")
                    nitrile_chemistry_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return nitrile_chemistry_detected
