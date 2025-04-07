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
    Detects if the synthesis route has an amide formation as one of the final steps (depth 0-1).
    """
    late_stage_amide_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_amide_formation

        if node["type"] == "reaction" and depth <= 1:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]
                reactants = reactants_smiles.split(".")

                # Check for amide formation reactions directly
                is_amide_formation = (
                    checker.check_reaction(
                        "Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi
                    )
                    or checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)
                    or checker.check_reaction(
                        "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
                    )
                    or checker.check_reaction("Acyl chloride with secondary amine to amide", rsmi)
                    or checker.check_reaction("Ester with primary amine to amide", rsmi)
                    or checker.check_reaction("Ester with secondary amine to amide", rsmi)
                    or checker.check_reaction("Acyl chloride with ammonia to amide", rsmi)
                    or checker.check_reaction("Schotten-Baumann_amide", rsmi)
                )

                # If direct reaction check fails, try checking reactants and products
                if not is_amide_formation:
                    # Check reactants for carboxylic acid/acyl chloride and amine
                    has_carboxylic_acid = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants
                    )
                    has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)
                    has_ester = any(checker.check_fg("Ester", r) for r in reactants)

                    has_primary_amine = any(checker.check_fg("Primary amine", r) for r in reactants)
                    has_secondary_amine = any(
                        checker.check_fg("Secondary amine", r) for r in reactants
                    )

                    # Check product for amide
                    has_amide_product = (
                        checker.check_fg("Primary amide", product_smiles)
                        or checker.check_fg("Secondary amide", product_smiles)
                        or checker.check_fg("Tertiary amide", product_smiles)
                    )

                    # Verify amide formation conditions
                    is_amide_formation = (
                        (has_carboxylic_acid or has_acyl_halide or has_ester)
                        and (has_primary_amine or has_secondary_amine)
                        and has_amide_product
                    )

                if is_amide_formation:
                    print(f"Found late-stage amide formation at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    late_stage_amide_formation = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return late_stage_amide_formation
