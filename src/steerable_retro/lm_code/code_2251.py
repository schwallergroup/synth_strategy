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
    Detects if the synthesis route involves a late-stage esterification
    (carboxylic acid to ester transformation in the last few steps).
    """
    esterification_found = False
    late_stage_threshold = 2  # Define what "late stage" means (depth < threshold)

    def dfs_traverse(node, depth=0):
        nonlocal esterification_found

        # Print current node for debugging
        if node["type"] == "mol":
            print(f"Examining molecule at depth {depth}: {node['smiles'][:30]}...")

        if node["type"] == "reaction" and depth < late_stage_threshold:
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi[:50]}...")

                # Check for various esterification reactions
                esterification_reactions = [
                    "Esterification of Carboxylic Acids",
                    "Schotten-Baumann to ester",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Oxidative esterification of primary alcohols",
                    "Transesterification",
                ]

                for rxn_type in esterification_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Found {rxn_type} reaction at depth {depth}")
                        esterification_found = True
                        return

                # Fallback check: look for carboxylic acid → ester transformation
                # In retrosynthesis, product contains acid and reactant contains ester
                product_has_acid = checker.check_fg("Carboxylic acid", product)
                reactant_has_ester = any(checker.check_fg("Ester", r) for r in reactants)

                # Also check the reverse (forward synthesis direction)
                reactant_has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)
                product_has_ester = checker.check_fg("Ester", product)

                if (product_has_acid and reactant_has_ester) or (
                    reactant_has_acid and product_has_ester
                ):
                    print(f"Found potential esterification at depth {depth} through FG analysis")

                    # Additional check to ensure it's not just coincidental presence of both groups
                    # Check if alcohol is also present in the reactants for forward direction
                    alcohol_types = [
                        "Primary alcohol",
                        "Secondary alcohol",
                        "Tertiary alcohol",
                        "Aromatic alcohol",
                    ]
                    reactant_has_alcohol = any(
                        any(checker.check_fg(alcohol, r) for alcohol in alcohol_types)
                        for r in reactants
                    )

                    if reactant_has_acid and product_has_ester and reactant_has_alcohol:
                        print(f"Confirmed esterification: acid + alcohol → ester at depth {depth}")
                        esterification_found = True
                        return

                    # For retrosynthetic direction
                    if product_has_acid and reactant_has_ester:
                        print(
                            f"Confirmed retrosynthetic esterification: ester → acid at depth {depth}"
                        )
                        esterification_found = True
                        return

            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Final result: esterification_found = {esterification_found}")
    return esterification_found
