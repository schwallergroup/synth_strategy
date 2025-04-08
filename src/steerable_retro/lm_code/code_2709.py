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
    Detects if the synthesis route involves Boc deprotection of a nitrogen.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check for any Boc deprotection reaction types
                boc_deprotection_types = [
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine",
                ]

                for rxn_type in boc_deprotection_types:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Detected {rxn_type} reaction")
                        result = True
                        return

                # If specific reaction check fails, try structural analysis
                # Check for Boc group in reactants but not in product
                has_boc_in_reactants = any(checker.check_fg("Boc", r) for r in reactants_smiles)
                has_boc_in_product = checker.check_fg("Boc", product_smiles)

                # Check for amine in both reactants and product
                has_amine_in_reactants = any(
                    checker.check_fg("Primary amine", r)
                    or checker.check_fg("Secondary amine", r)
                    or checker.check_fg("Tertiary amine", r)
                    for r in reactants_smiles
                )
                has_amine_in_product = (
                    checker.check_fg("Primary amine", product_smiles)
                    or checker.check_fg("Secondary amine", product_smiles)
                    or checker.check_fg("Tertiary amine", product_smiles)
                )

                if (
                    has_boc_in_reactants
                    and not has_boc_in_product
                    and has_amine_in_reactants
                    and has_amine_in_product
                ):
                    print("Detected Boc deprotection through structural analysis")
                    result = True
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal if result is still False
        if not result:
            for child in node.get("children", []):
                dfs_traverse(child)

    dfs_traverse(route)
    return result
