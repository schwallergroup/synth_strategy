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
    Checks if the route contains Boc protection and deprotection steps.
    """
    has_protection = False
    has_deprotection = False

    def dfs(node, depth=0):
        nonlocal has_protection, has_deprotection

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]

            # Check for Boc protection reactions
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction(
                    "Boc amine protection with Boc anhydride", rsmi
                )
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction(
                    "Boc amine protection of secondary amine", rsmi
                )
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                has_protection = True
                print(f"Found Boc protection reaction: {rsmi}")

            # Check for Boc deprotection reactions
            if (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            ):
                has_deprotection = True
                print(f"Found Boc deprotection reaction: {rsmi}")

            # If both protection and deprotection are found, return True
            if has_protection and has_deprotection:
                return True

        # Check children
        for child in node.get("children", []):
            if dfs(child, depth + 1):
                return True

        return has_protection and has_deprotection

    return dfs(route)
