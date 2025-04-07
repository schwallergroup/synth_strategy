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
    This function detects if the synthesis uses a late-stage amide formation
    as the final or penultimate step.
    """
    amide_formation_depths = []

    def dfs_traverse(node):
        if node["type"] == "reaction":
            depth = node.get("metadata", {}).get("depth", 0)
            rsmi = node.get("metadata", {}).get("rsmi", "")

            if not rsmi:
                return

            # Check for various amide formation reactions
            amide_formation_reactions = [
                "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                "Carboxylic acid with primary amine to amide",
                "Acyl chloride with ammonia to amide",
                "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                "Acyl chloride with secondary amine to amide",
                "Ester with ammonia to amide",
                "Ester with primary amine to amide",
                "Ester with secondary amine to amide",
                "Schotten-Baumann_amide",
            ]

            for reaction_type in amide_formation_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found amide formation reaction: {reaction_type} at depth {depth}")
                    amide_formation_depths.append(depth)
                    break

            # If no specific reaction type matched, check for amide formation by functional group changes
            if not any(checker.check_reaction(r, rsmi) for r in amide_formation_reactions):
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants contain carboxylic acid or acyl halide and amine
                has_carboxylic_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                has_ester = any(checker.check_fg("Ester", r) for r in reactants)
                has_amine = any(checker.check_fg("Primary amine", r) for r in reactants) or any(
                    checker.check_fg("Secondary amine", r) for r in reactants
                )

                # Check if product contains amide
                has_amide_product = (
                    checker.check_fg("Primary amide", product)
                    or checker.check_fg("Secondary amide", product)
                    or checker.check_fg("Tertiary amide", product)
                )

                if (
                    (has_carboxylic_acid or has_acyl_halide or has_ester)
                    and has_amine
                    and has_amide_product
                ):
                    print(f"Found amide formation by functional group analysis at depth {depth}")
                    amide_formation_depths.append(depth)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Check if amide formation occurs at depth 0 or 1 (late stage)
    if any(depth <= 1 for depth in amide_formation_depths):
        print(
            f"Found late-stage amide formation at depths: {[d for d in amide_formation_depths if d <= 1]}"
        )
        return True
    return False
