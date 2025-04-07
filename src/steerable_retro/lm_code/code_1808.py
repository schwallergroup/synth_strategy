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
    Detects if the synthesis route involves late-stage amide formation.
    Late-stage means it occurs at a low depth in the synthesis tree.
    """
    amide_formation_depths = []

    def dfs(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]

                # Check for amide formation reactions
                if (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                        rsmi,
                    )
                    or checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with ammonia to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Carboxylic acid to amide conversion", rsmi)
                    or checker.check_reaction("Schotten-Baumann to ester", rsmi)
                    or checker.check_reaction("{Schotten-Baumann_amide}", rsmi)
                ):

                    # Verify that this reaction actually forms an amide
                    product = rsmi.split(">")[-1]

                    # Check if product has an amide group
                    has_amide_product = (
                        checker.check_fg("Primary amide", product)
                        or checker.check_fg("Secondary amide", product)
                        or checker.check_fg("Tertiary amide", product)
                    )

                    if has_amide_product:
                        amide_formation_depths.append(depth)
                        print(f"Found amide formation at depth {depth}: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction node for amide formation: {e}")

        # Recursively process children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS traversal
    dfs(route)

    # Check if we found any amide formation reactions
    if not amide_formation_depths:
        print("No amide formation reactions found")
        return False

    # Check if at least one amide formation occurs at a late stage (low depth)
    # We'll consider depth <= 2 as late stage
    late_stage = any(depth <= 2 for depth in amide_formation_depths)

    print(f"Late-stage amide formation: {late_stage}")
    return late_stage
