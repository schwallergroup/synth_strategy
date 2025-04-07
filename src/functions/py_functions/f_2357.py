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
    Detects silyl protection-deprotection strategy in the synthesis route.
    """
    protection_found = False
    deprotection_found = False

    def dfs(node, depth=0):
        nonlocal protection_found, deprotection_found

        if node["type"] == "reaction":
            try:
                rxn_smiles = node.get("metadata", {}).get("rsmi", "")
                if not rxn_smiles:
                    return

                # Check for silyl protection reaction
                if checker.check_reaction(
                    "Alcohol protection with silyl ethers", rxn_smiles
                ):
                    protection_found = True
                    print(f"Found silyl protection reaction: {rxn_smiles}")

                # Check for silyl deprotection reactions
                if (
                    checker.check_reaction(
                        "Alcohol deprotection from silyl ethers", rxn_smiles
                    )
                    or checker.check_reaction(
                        "Alcohol deprotection from silyl ethers (double)", rxn_smiles
                    )
                    or checker.check_reaction(
                        "Alcohol deprotection from silyl ethers (diol)", rxn_smiles
                    )
                ):
                    deprotection_found = True
                    print(f"Found silyl deprotection reaction: {rxn_smiles}")

                # Additional check for TMS protection/deprotection
                reactants = rxn_smiles.split(">")[0].split(".")
                product = rxn_smiles.split(">")[-1]

                # Check for TMS ether protective group in reactants or products
                if any(
                    checker.check_fg("TMS ether protective group", r) for r in reactants
                ) or checker.check_fg("TMS ether protective group", product):
                    if any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants
                    ):
                        protection_found = True
                        print(f"Found silyl protection via TMS group: {rxn_smiles}")
                    elif (
                        checker.check_fg("Primary alcohol", product)
                        or checker.check_fg("Secondary alcohol", product)
                        or checker.check_fg("Tertiary alcohol", product)
                    ):
                        deprotection_found = True
                        print(f"Found silyl deprotection to alcohol: {rxn_smiles}")

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    print(
        f"Silyl protection found: {protection_found}, deprotection found: {deprotection_found}"
    )
    # Consider the strategy detected if either protection or deprotection is found
    return protection_found or deprotection_found
