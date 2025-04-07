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
    This function detects Williamson ether synthesis (phenol + alkyl halide → ether).
    """
    williamson_detected = False

    def dfs_traverse(node):
        nonlocal williamson_detected

        if williamson_detected:
            return  # Early return if already detected

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]

            # First, check if this is a Williamson ether synthesis reaction using predefined patterns
            if checker.check_reaction("Williamson Ether Synthesis", rsmi):
                print("Detected Williamson ether synthesis via reaction pattern")
                williamson_detected = True
                return

            # Also check for intramolecular variant
            if checker.check_reaction("Williamson Ether Synthesis (intra to epoxy)", rsmi):
                print("Detected intramolecular Williamson ether synthesis via reaction pattern")
                williamson_detected = True
                return

            # If the direct check fails, verify the reaction components
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Analyzing reaction: {rsmi}")

            # Check for phenol in reactants
            has_phenol = any(checker.check_fg("Phenol", r) for r in reactants)

            # Check for various types of halides in reactants
            has_primary_halide = any(checker.check_fg("Primary halide", r) for r in reactants)
            has_secondary_halide = any(checker.check_fg("Secondary halide", r) for r in reactants)
            has_tertiary_halide = any(checker.check_fg("Tertiary halide", r) for r in reactants)
            has_aromatic_halide = any(checker.check_fg("Aromatic halide", r) for r in reactants)

            has_alkyl_halide = has_primary_halide or has_secondary_halide or has_tertiary_halide
            has_halide = has_alkyl_halide or has_aromatic_halide

            # Check for ether in product
            has_ether = checker.check_fg("Ether", product)

            # Also check for alcohol in reactants (for alternative mechanism)
            has_alcohol = (
                any(checker.check_fg("Primary alcohol", r) for r in reactants)
                or any(checker.check_fg("Secondary alcohol", r) for r in reactants)
                or any(checker.check_fg("Tertiary alcohol", r) for r in reactants)
                or any(checker.check_fg("Aromatic alcohol", r) for r in reactants)
            )

            print(
                f"Components - Phenol: {has_phenol}, Alkyl halide: {has_alkyl_halide}, "
                f"Aromatic halide: {has_aromatic_halide}, Ether in product: {has_ether}, "
                f"Alcohol: {has_alcohol}"
            )

            # Verify this is a Williamson ether synthesis
            if has_phenol and has_halide and has_ether:
                print(f"Detected Williamson ether synthesis: Phenol + Halide → Ether")
                williamson_detected = True
                return

            # Check for alternative mechanism: alcohol + halide → ether
            if has_alcohol and has_halide and has_ether:
                print(f"Detected Williamson-like ether synthesis: Alcohol + Halide → Ether")
                williamson_detected = True
                return

            # Check for Mitsunobu reaction which can also form ethers
            if checker.check_reaction("Mitsunobu aryl ether", rsmi):
                print("Detected Mitsunobu aryl ether (similar to Williamson)")
                williamson_detected = True
                return

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if not williamson_detected:
        print("No Williamson ether synthesis detected in the route")

    return williamson_detected
