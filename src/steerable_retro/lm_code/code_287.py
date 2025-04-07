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
    This function detects if the synthesis route uses late-stage amide formation
    as the final step from a carboxylic acid.
    """
    late_stage_amide = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide

        print(f"Examining node at depth {depth}, type: {node['type']}")

        # Check if this is a reaction node at depth 0 or 1 (late stage)
        if node["type"] == "reaction" and depth <= 1:
            print(f"Examining potential late-stage reaction at depth {depth}")

            if "metadata" in node and "rsmi" in node["metadata"]:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    print(f"Reaction SMILES: {rsmi}")

                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Reactants: {reactants}")
                    print(f"Product: {product}")

                    # Check if this is an amide formation reaction directly
                    amide_reaction_types = [
                        "Carboxylic acid with primary amine to amide",
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                        "Acyl chloride with secondary amine to amide",
                        "Ester with primary amine to amide",
                        "Ester with secondary amine to amide",
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                        "Acylation of primary amines",
                        "Acylation of secondary amines",
                        "Schotten-Baumann_amide",
                        "Acyl chloride with ammonia to amide",
                    ]

                    is_amide_formation = False
                    for reaction_type in amide_reaction_types:
                        if checker.check_reaction(reaction_type, rsmi):
                            is_amide_formation = True
                            print(f"Detected amide formation reaction type: {reaction_type}")
                            break

                    # If direct reaction check failed, check for functional groups
                    if not is_amide_formation:
                        # Check if reactants contain carboxylic acid or acyl halide
                        has_acid = False
                        has_acyl_halide = False
                        has_amine = False

                        for reactant in reactants:
                            if checker.check_fg("Carboxylic acid", reactant):
                                has_acid = True
                                print(f"Found carboxylic acid in reactant: {reactant}")

                            if checker.check_fg("Acyl halide", reactant):
                                has_acyl_halide = True
                                print(f"Found acyl halide in reactant: {reactant}")

                            if (
                                checker.check_fg("Primary amine", reactant)
                                or checker.check_fg("Secondary amine", reactant)
                                or checker.check_fg("Aniline", reactant)
                            ):
                                has_amine = True
                                print(f"Found amine in reactant: {reactant}")

                        # Check if product contains amide
                        has_amide = False
                        if (
                            checker.check_fg("Primary amide", product)
                            or checker.check_fg("Secondary amide", product)
                            or checker.check_fg("Tertiary amide", product)
                        ):
                            has_amide = True
                            print(f"Found amide in product: {product}")

                        # If we have the right functional groups, check if this could be amide formation
                        if has_amide and (
                            (has_acid and has_amine) or (has_acyl_halide and has_amine)
                        ):
                            is_amide_formation = True
                            print("Detected amide formation based on functional groups")

                    if is_amide_formation:
                        late_stage_amide = True
                        print(f"Detected late-stage amide formation at depth {depth}")

                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: {late_stage_amide}")
    return late_stage_amide
