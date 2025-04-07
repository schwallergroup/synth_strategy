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
    This function detects if the synthesis has a late-stage amide coupling
    (amide formation in the final step or penultimate step).
    """
    has_late_stage_amide = False

    # Determine if route starts with a molecule (product) or reaction
    starts_with_molecule = route["type"] == "mol"

    def is_amide_formation_reaction(rsmi, reactants_smiles, product_smiles):
        """Helper function to check if a reaction is an amide formation"""
        # Check for amide in product
        has_amide_in_product = (
            checker.check_fg("Primary amide", product_smiles)
            or checker.check_fg("Secondary amide", product_smiles)
            or checker.check_fg("Tertiary amide", product_smiles)
        )

        print(f"Product contains amide: {has_amide_in_product}")

        if not has_amide_in_product:
            return False

        # Check for amide in reactants
        has_amide_in_reactants = any(
            checker.check_fg("Primary amide", r)
            or checker.check_fg("Secondary amide", r)
            or checker.check_fg("Tertiary amide", r)
            for r in reactants_smiles
        )

        print(f"Reactants contain amide: {has_amide_in_reactants}")

        if has_amide_in_reactants:
            return False  # Amide was already present, not newly formed

        # Check for amide coupling reaction types
        is_amide_coupling = (
            checker.check_reaction(
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
            )
            or checker.check_reaction(
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                rsmi,
            )
            or checker.check_reaction(
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                rsmi,
            )
            or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
            or checker.check_reaction(
                "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
            )
            or checker.check_reaction(
                "Carboxylic acid with primary amine to amide", rsmi
            )
            or checker.check_reaction("Ester with primary amine to amide", rsmi)
            or checker.check_reaction("Ester with ammonia to amide", rsmi)
            or checker.check_reaction(
                "Acyl chloride with secondary amine to amide", rsmi
            )
            or checker.check_reaction("Ester with secondary amine to amide", rsmi)
            or checker.check_reaction("Acylation of primary amines", rsmi)
            or checker.check_reaction("Acylation of secondary amines", rsmi)
            or checker.check_reaction("Schotten-Baumann_amide", rsmi)
            or checker.check_reaction("Schotten-Baumann to ester", rsmi)
            or checker.check_reaction("Carboxylic acid to amide conversion", rsmi)
        )

        print(f"Reaction is amide coupling: {is_amide_coupling}")

        if is_amide_coupling:
            return True

        # Alternative check using reactant functional groups if reaction check fails
        has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants_smiles)
        has_acyl_halide = any(
            checker.check_fg("Acyl halide", r) for r in reactants_smiles
        )
        has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)
        has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants_smiles)

        has_primary_amine = any(
            checker.check_fg("Primary amine", r) for r in reactants_smiles
        )
        has_secondary_amine = any(
            checker.check_fg("Secondary amine", r) for r in reactants_smiles
        )
        has_ammonia = any(
            "N" in r and len(r) <= 3 for r in reactants_smiles
        )  # Simple check for ammonia

        if (has_acid or has_acyl_halide or has_ester or has_anhydride) and (
            has_primary_amine or has_secondary_amine or has_ammonia
        ):
            print("Detected amide coupling through functional group analysis")
            return True

        return False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_stage_amide

        print(f"Traversing node type: {node['type']} at depth: {depth}")

        # Check if this is a late-stage reaction
        is_late_stage = (starts_with_molecule and depth <= 2) or (
            not starts_with_molecule and depth <= 1
        )

        if node["type"] == "reaction" and is_late_stage:
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Analyzing potential late-stage reaction: {rsmi}")

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                if is_amide_formation_reaction(rsmi, reactants_smiles, product_smiles):
                    has_late_stage_amide = True
                    print("Detected late-stage amide coupling")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_late_stage_amide
