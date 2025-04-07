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
    This function detects a synthesis strategy involving multiple C-O bond formations.
    """
    co_formation_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Depth {depth}, Examining reaction: {rsmi}")

            # Check for various C-O bond formation reactions
            co_formation_reaction_types = [
                "Williamson Ether Synthesis",
                "Esterification of Carboxylic Acids",
                "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
                "Schotten-Baumann to ester",
                "Transesterification",
                "Alcohol protection with silyl ethers",
                "O-alkylation of carboxylic acids with diazo compounds",
                "Formation of Sulfonic Esters",
                "Oxidative esterification of primary alcohols",
                "Mitsunobu esterification",
                "Mitsunobu aryl ether",
                "Chan-Lam alcohol",
                "Chan-Lam etherification",
                "Acetic anhydride and alcohol to ester",
                "Alcohol to ether",
                "Williamson Ether Synthesis (intra to epoxy)",
                "Oxidation of alcohol to carboxylic acid",
                "Oxidation of aldehydes to carboxylic acids",
                "Oxidation of alcohol and aldehyde to ester",
            ]

            found_co_formation = False
            for reaction_type in co_formation_reaction_types:
                if checker.check_reaction(reaction_type, rsmi):
                    print(f"Found C-O bond formation: {reaction_type}")
                    co_formation_reactions.append((depth, reaction_type, rsmi))
                    found_co_formation = True
                    break

            # If no specific reaction type was found, check for C-O bond formation by examining
            # functional groups in reactants and products
            if not found_co_formation:
                try:
                    reactants_part = rsmi.split(">")[0]
                    products_part = rsmi.split(">")[-1]

                    # Check for alcohol, ether, or ester formation
                    reactants_have_alcohol = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants_part.split(".")
                    )

                    product_has_ether = checker.check_fg("Ether", products_part)
                    product_has_ester = checker.check_fg("Ester", products_part)

                    if reactants_have_alcohol and (product_has_ether or product_has_ester):
                        print(f"Found C-O bond formation through FG analysis")
                        co_formation_reactions.append((depth, "C-O bond formation", rsmi))
                except Exception as e:
                    print(f"Error in FG analysis: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Total C-O bond formations found: {len(co_formation_reactions)}")
    for depth, reaction_type, rsmi in co_formation_reactions:
        print(f"Depth {depth}: {reaction_type} - {rsmi}")

    # Check if we have at least 2 C-O bond formations
    if len(co_formation_reactions) < 2:
        return False

    # Sort reactions by depth to check if they're sequential
    co_formation_reactions.sort(key=lambda x: x[0])

    # Check if the reactions are in different branches or in sequence
    # For simplicity, we'll consider them sequential if they occur at different depths
    depths = [d for d, _, _ in co_formation_reactions]
    unique_depths = set(depths)

    print(f"Reaction depths: {depths}")
    print(f"Unique depths: {unique_depths}")

    # If we have at least 2 different depths, consider it sequential
    return len(unique_depths) >= 2
