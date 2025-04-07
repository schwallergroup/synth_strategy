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
    Detects if the synthesis route involves reductive amination reactions.
    """
    reductive_amination_reactions = []

    def find_reductive_amination(node, depth=0):
        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for reductive amination reactions
            if (
                checker.check_reaction("Reductive amination with aldehyde", rxn_smiles)
                or checker.check_reaction("Reductive amination with ketone", rxn_smiles)
                or checker.check_reaction(
                    "Reductive amination with alcohol", rxn_smiles
                )
                or checker.check_reaction("reductive amination", rxn_smiles)
                or checker.check_reaction("{reductive amination}", rxn_smiles)
                or checker.check_reaction("Mignonac reaction", rxn_smiles)
            ):
                reductive_amination_reactions.append((rxn_smiles, depth))

            # Check for functional group transformations that might indicate reductive amination
            reactants = rxn_smiles.split(">")[0].split(".")
            product = rxn_smiles.split(">")[-1]

            # Check if reactants contain aldehyde/ketone and amine, and product contains amine
            has_carbonyl = any(
                checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                for r in reactants
            )
            has_amine = any(
                checker.check_fg("Primary amine", r)
                or checker.check_fg("Secondary amine", r)
                for r in reactants
            )
            product_has_amine = checker.check_fg(
                "Secondary amine", product
            ) or checker.check_fg("Tertiary amine", product)

            if has_carbonyl and has_amine and product_has_amine:
                reductive_amination_reactions.append((rxn_smiles, depth))

        for child in node.get("children", []):
            find_reductive_amination(child, depth + 1)

    find_reductive_amination(route)

    # Remove duplicates
    unique_reactions = set(rxn for rxn, _ in reductive_amination_reactions)

    # Consider it a reductive amination strategy if at least one reductive amination is found
    print(f"Reductive amination reactions found: {len(unique_reactions)}")
    return len(unique_reactions) > 0
