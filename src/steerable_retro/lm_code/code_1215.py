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
    """Check for Boc protection followed by deprotection in the synthesis route"""
    has_protection = False
    has_deprotection = False

    def dfs(node, depth=0):
        nonlocal has_protection, has_deprotection

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for Boc protection reaction
            if (
                checker.check_reaction("Boc amine protection", rxn_smiles)
                or checker.check_reaction("Boc amine protection explicit", rxn_smiles)
                or checker.check_reaction("Boc amine protection with Boc anhydride", rxn_smiles)
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rxn_smiles)
                or checker.check_reaction("Boc amine protection of secondary amine", rxn_smiles)
                or checker.check_reaction("Boc amine protection of primary amine", rxn_smiles)
            ):
                has_protection = True
                print(f"Found Boc protection at depth {depth}: {rxn_smiles}")

            # Check for Boc deprotection reaction
            if (
                checker.check_reaction("Boc amine deprotection", rxn_smiles)
                or checker.check_reaction("Boc amine deprotection of guanidine", rxn_smiles)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rxn_smiles)
                or checker.check_reaction("Tert-butyl deprotection of amine", rxn_smiles)
            ):
                has_deprotection = True
                print(f"Found Boc deprotection at depth {depth}: {rxn_smiles}")

            # Additional check for Boc-containing molecules in reactants
            if not has_protection:
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    for r in reactants:
                        if "BOC" in r.upper() or "OC(C)(C)C" in r:
                            has_protection = True
                            print(f"Found Boc-containing reactant at depth {depth}: {r}")
                            break
                except Exception as e:
                    print(f"Error checking Boc-containing reactants: {e}")

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return has_protection or has_deprotection  # Changed to OR instead of AND to be less restrictive
