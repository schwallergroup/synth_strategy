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
    This function detects if the synthetic route involves late-stage amide formation (at depth 0 or 1).
    """
    late_stage_amide = False

    def dfs_traverse(node, current_depth=0):
        nonlocal late_stage_amide

        if node["type"] == "reaction" and current_depth <= 1:  # Late stage (depth 0 or 1)
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for amide formation reaction types
                amide_formation_reactions = [
                    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
                    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
                    "Acyl chloride with ammonia to amide",
                    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
                    "Acyl chloride with secondary amine to amide",
                    "Carboxylic acid with primary amine to amide",
                    "Ester with ammonia to amide",
                    "Ester with primary amine to amide",
                    "Ester with secondary amine to amide",
                    "Schotten-Baumann_amide",
                    "Acylation of primary amines",
                    "Acylation of secondary amines",
                    "Carboxylic acid to amide conversion",
                ]

                is_amide_formation = any(
                    checker.check_reaction(rxn, rsmi) for rxn in amide_formation_reactions
                )

                # If reaction type check fails, check for amide formation by structure
                if not is_amide_formation:
                    # Count amides in reactants and product
                    amide_types = ["Primary amide", "Secondary amide", "Tertiary amide"]

                    # Check product for amides
                    product_amide_count = sum(
                        1
                        for amide_type in amide_types
                        if checker.check_fg(amide_type, product_smiles)
                    )

                    # Check reactants for amides
                    reactants_amide_count = 0
                    for reactant in reactants_smiles:
                        reactants_amide_count += sum(
                            1
                            for amide_type in amide_types
                            if checker.check_fg(amide_type, reactant)
                        )

                    # If product has more amides than reactants, amide formation occurred
                    is_amide_formation = product_amide_count > reactants_amide_count

                if is_amide_formation:
                    print(f"Found late-stage amide formation at depth {current_depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    late_stage_amide = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage amide formation strategy detected: {late_stage_amide}")
    return late_stage_amide
