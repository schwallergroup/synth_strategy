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
    This function detects if the synthetic route involves Boc protection of an amine.
    """
    boc_protection_found = False

    def dfs_traverse(node):
        nonlocal boc_protection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a Boc protection reaction using the checker function
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rsmi)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
                or checker.check_reaction("Boc amine protection of secondary amine", rsmi)
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
            ):
                print(f"Boc protection reaction detected: {rsmi}")
                boc_protection_found = True
            else:
                # Fallback check by examining functional groups
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if reactants contain an amine and product contains a Boc group
                has_primary_amine = any(
                    checker.check_fg("Primary amine", r) for r in reactants if r
                )
                has_secondary_amine = any(
                    checker.check_fg("Secondary amine", r) for r in reactants if r
                )
                has_boc_protected_amine = checker.check_fg("Boc", product) if product else False

                if (has_primary_amine or has_secondary_amine) and has_boc_protected_amine:
                    print(f"Boc protection detected by functional group analysis: {rsmi}")
                    boc_protection_found = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_protection_found
