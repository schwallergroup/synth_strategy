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
    This function detects if the synthetic route involves an esterification step.
    """
    has_esterification = False

    def dfs_traverse(node):
        nonlocal has_esterification

        if node["type"] == "reaction" and not has_esterification:
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Direct reaction type checking - most reliable method
                esterification_reactions = [
                    "Esterification of Carboxylic Acids",
                    "Transesterification",
                    "Schotten-Baumann to ester",
                    "O-alkylation of carboxylic acids with diazo compounds",
                    "Oxidative esterification of primary alcohols",
                ]

                for rxn_type in esterification_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        has_esterification = True
                        print(f"Found esterification reaction: {rxn_type}")
                        break

                # Fallback to functional group analysis if reaction type check fails
                if not has_esterification:
                    # Check for carboxylic acid in reactants
                    has_acid = any(checker.check_fg("Carboxylic acid", r) for r in reactants)

                    # Check for alcohol in reactants
                    has_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants
                    )

                    # Check for ester in product
                    has_ester_product = checker.check_fg("Ester", product)

                    # Check for acyl halides which can also form esters
                    has_acyl_halide = any(checker.check_fg("Acyl halide", r) for r in reactants)

                    # Check for anhydrides which can also form esters
                    has_anhydride = any(checker.check_fg("Anhydride", r) for r in reactants)

                    if has_ester_product and (
                        (has_acid and has_alcohol)
                        or (has_acyl_halide and has_alcohol)
                        or (has_anhydride and has_alcohol)
                    ):
                        has_esterification = True
                        print("Found esterification by functional group analysis")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_esterification
