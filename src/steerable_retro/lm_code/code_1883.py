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
    Detects if the route contains a Grignard reaction (organometallic addition to carbonyl)
    """
    has_grignard = False

    def dfs_traverse(node):
        nonlocal has_grignard

        if node["type"] == "reaction" and "metadata" in node:
            # Safely extract reaction SMILES
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            # Check for various Grignard reaction types using the checker function
            grignard_reactions = [
                "Grignard from aldehyde to alcohol",
                "Grignard from ketone to alcohol",
                "Formation of Grignard reagents",
                "Grignard with CO2 to carboxylic acid",
                "Olefination of ketones with Grignard reagents",
                "Olefination of aldehydes with Grignard reagents",
                "Grignard from nitrile to ketone",
                "Preparation of trialkylsilanes with Grignard reagents",
                "Grignard_carbonyl",
                "Grignard_alcohol",
            ]

            for reaction_type in grignard_reactions:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found Grignard reaction: {reaction_type}")
                    print(f"Reaction SMILES: {rsmi}")
                    has_grignard = True
                    return  # Exit early once found

            # If no direct reaction match, check for characteristic functional groups
            try:
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant contains a magnesium halide (characteristic of Grignard reagents)
                has_mg_halide = any(
                    checker.check_fg("Magnesium halide", reactant) for reactant in reactants
                )

                # Check if any reactant contains a carbonyl group that Grignard reagents typically react with
                has_carbonyl = any(
                    checker.check_fg("Aldehyde", reactant)
                    or checker.check_fg("Ketone", reactant)
                    or checker.check_fg("Ester", reactant)
                    for reactant in reactants
                )

                # If we have both a Grignard reagent and a suitable carbonyl, it's likely a Grignard reaction
                if has_mg_halide and has_carbonyl:
                    print(f"Found Grignard reaction based on functional groups")
                    print(f"Reaction SMILES: {rsmi}")
                    has_grignard = True
                    return
            except Exception as e:
                print(f"Error analyzing reaction components: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_grignard
