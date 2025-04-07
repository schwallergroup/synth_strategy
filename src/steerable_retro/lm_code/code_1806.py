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
    Detects if the synthesis route involves Boc protection followed by deprotection.
    """
    # Track protection and deprotection reactions
    protection_reactions = []
    deprotection_reactions = []

    def dfs(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]

                # Check for Boc protection reaction
                if (
                    checker.check_reaction("Boc amine protection", rsmi)
                    or checker.check_reaction("Boc amine protection explicit", rsmi)
                    or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                    or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                    or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                    or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                ):
                    protection_reactions.append((depth, rsmi))
                    print(f"Found Boc protection reaction at depth {depth}: {rsmi}")

                # Check for Boc deprotection reaction
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    deprotection_reactions.append((depth, rsmi))
                    print(f"Found Boc deprotection reaction at depth {depth}: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction node for Boc sequence: {e}")

        # Recursively process children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    # Start DFS traversal
    dfs(route)

    # Check if we found both protection and deprotection
    has_protection = len(protection_reactions) > 0
    has_deprotection = len(deprotection_reactions) > 0

    # Check if protection happens before deprotection (higher depth)
    correct_sequence = False
    if has_protection and has_deprotection:
        # Get the minimum depth for each (earliest occurrence)
        min_protection_depth = min([d for d, _ in protection_reactions])
        min_deprotection_depth = min([d for d, _ in deprotection_reactions])

        # Protection should happen at a higher depth (earlier in synthesis)
        correct_sequence = min_protection_depth > min_deprotection_depth

    result = has_protection and has_deprotection and correct_sequence
    print(f"Boc protection/deprotection sequence: {result}")
    return result
