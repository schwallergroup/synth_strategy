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
    Checks if the synthetic route contains ester interconversion reactions.
    """
    ester_formation_reactions = []
    ester_hydrolysis_reactions = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for ester formation reactions
            if (
                checker.check_reaction("Esterification of Carboxylic Acids", rxn_smiles)
                or checker.check_reaction("Schotten-Baumann to ester", rxn_smiles)
                or checker.check_reaction(
                    "O-alkylation of carboxylic acids with diazo compounds", rxn_smiles
                )
                or checker.check_reaction(
                    "Oxidative esterification of primary alcohols", rxn_smiles
                )
                or checker.check_reaction("Transesterification", rxn_smiles)
                or checker.check_reaction("Acetic anhydride and alcohol to ester", rxn_smiles)
                or checker.check_reaction("Mitsunobu esterification", rxn_smiles)
                or checker.check_reaction(
                    "Intramolecular transesterification/Lactone formation", rxn_smiles
                )
            ):
                ester_formation_reactions.append((depth, rxn_smiles))
                print(f"Ester formation reaction found at depth {depth}: {rxn_smiles[:50]}...")

            # Check for ester hydrolysis reactions
            if (
                checker.check_reaction(
                    "Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rxn_smiles
                )
                or checker.check_reaction("Ester saponification (methyl deprotection)", rxn_smiles)
                or checker.check_reaction("Ester saponification (alkyl deprotection)", rxn_smiles)
                or checker.check_reaction("COOH ethyl deprotection", rxn_smiles)
            ):
                ester_hydrolysis_reactions.append((depth, rxn_smiles))
                print(f"Ester hydrolysis reaction found at depth {depth}: {rxn_smiles[:50]}...")

            # Also check for functional group changes that might indicate ester interconversion
            reactants = rxn_smiles.split(">")[0].split(".")
            product = rxn_smiles.split(">")[-1]

            # Check for ester formation patterns
            if any(checker.check_fg("Carboxylic acid", r) for r in reactants) and checker.check_fg(
                "Ester", product
            ):
                ester_formation_reactions.append((depth, rxn_smiles))
                print(f"Carboxylic acid to ester pattern found at depth {depth}")

            # Check for ester hydrolysis patterns
            if any(checker.check_fg("Ester", r) for r in reactants) and checker.check_fg(
                "Carboxylic acid", product
            ):
                ester_hydrolysis_reactions.append((depth, rxn_smiles))
                print(f"Ester to carboxylic acid pattern found at depth {depth}")

        # Recursively traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both ester formation and hydrolysis reactions
    if ester_formation_reactions and ester_hydrolysis_reactions:
        # Check if any formation happens before any hydrolysis
        valid_sequence = any(
            f_depth < h_depth
            for f_depth, _ in ester_formation_reactions
            for h_depth, _ in ester_hydrolysis_reactions
        )

        if valid_sequence:
            print(
                f"Ester interconversion detected: {len(ester_formation_reactions)} formation reactions, {len(ester_hydrolysis_reactions)} hydrolysis reactions"
            )
            return True
        else:
            print(f"Ester reactions not in correct order")
    else:
        print(
            f"Ester interconversion not complete. Formation reactions: {len(ester_formation_reactions)}, Hydrolysis reactions: {len(ester_hydrolysis_reactions)}"
        )

    # If we have at least one ester formation or hydrolysis, consider it partial evidence
    if ester_formation_reactions or ester_hydrolysis_reactions:
        print("Partial ester interconversion detected")
        return True

    return False
