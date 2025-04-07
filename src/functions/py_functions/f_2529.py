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
    This function detects a strategy involving halogenation (fluorination, chlorination,
    bromination, or iodination) as a key step in the synthesis.
    """
    # Track if we found halogenation
    halogenation_found = False

    def dfs_traverse(node):
        nonlocal halogenation_found

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]

                try:
                    # Check for halogenation reactions using the checker function
                    if (
                        checker.check_reaction("Aromatic fluorination", rsmi)
                        or checker.check_reaction("Aromatic chlorination", rsmi)
                        or checker.check_reaction("Aromatic bromination", rsmi)
                        or checker.check_reaction("Aromatic iodination", rsmi)
                        or checker.check_reaction("Chlorination", rsmi)
                        or checker.check_reaction("Fluorination", rsmi)
                        or checker.check_reaction("Bromination", rsmi)
                        or checker.check_reaction("Iodination", rsmi)
                    ):

                        print(f"Detected halogenation reaction: {rsmi}")
                        halogenation_found = True

                    # If no direct halogenation reaction is detected, check for halogen addition
                    if not halogenation_found:
                        reactants_smiles = rsmi.split(">")[0].split(".")
                        product_smiles = rsmi.split(">")[-1]

                        # Check if product contains halogen functional groups that weren't in reactants
                        product_has_halogen = (
                            checker.check_fg("Aromatic halide", product_smiles)
                            or checker.check_fg("Primary halide", product_smiles)
                            or checker.check_fg("Secondary halide", product_smiles)
                            or checker.check_fg("Tertiary halide", product_smiles)
                            or checker.check_fg("Alkenyl halide", product_smiles)
                            or checker.check_fg("Haloalkyne", product_smiles)
                        )

                        if product_has_halogen:
                            # Check if any reactant has the same halogen functional group
                            reactants_have_halogen = any(
                                checker.check_fg("Aromatic halide", r)
                                or checker.check_fg("Primary halide", r)
                                or checker.check_fg("Secondary halide", r)
                                or checker.check_fg("Tertiary halide", r)
                                or checker.check_fg("Alkenyl halide", r)
                                or checker.check_fg("Haloalkyne", r)
                                for r in reactants_smiles
                            )

                            # If product has halogen but reactants don't, it's a halogenation
                            if not reactants_have_halogen:
                                print(f"Detected halogen addition: {rsmi}")
                                halogenation_found = True

                except Exception as e:
                    print(f"Error processing reaction SMILES: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return halogenation_found
