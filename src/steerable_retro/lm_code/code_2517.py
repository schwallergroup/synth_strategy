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
    This function detects if the synthetic route involves Boc protection/deprotection of an amine.
    """
    boc_protected_found = False
    protection_found = False
    deprotection_found = False

    def dfs_traverse(node):
        nonlocal boc_protected_found, protection_found, deprotection_found

        if node["type"] == "mol":
            # Check if molecule contains Boc-protected amine
            if node.get("smiles"):
                mol_smiles = node["smiles"]
                # Exclude Boc anhydride which is not a Boc-protected amine
                if (
                    checker.check_fg("Boc", mol_smiles)
                    and mol_smiles != "CC(C)(C)OC(=O)OC(=O)OC(C)(C)C"
                ):
                    boc_protected_found = True
                    print(f"Boc-protected amine found in molecule: {mol_smiles}")

        elif node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
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
                    protection_found = True
                    print(f"Boc protection reaction detected: {rsmi}")

                # Check for Boc deprotection reaction
                if (
                    checker.check_reaction("Boc amine deprotection", rsmi)
                    or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                    or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                    or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
                ):
                    deprotection_found = True
                    print(f"Boc deprotection reaction detected: {rsmi}")

                # If reaction type checkers didn't find anything, check manually
                if not protection_found and not deprotection_found:
                    reactants = rsmi.split(">")[0].split(".")
                    products = rsmi.split(">")[-1].split(".")

                    # Check if Boc group is added (protection)
                    reactants_have_boc = any(checker.check_fg("Boc", r) for r in reactants)
                    products_have_boc = any(checker.check_fg("Boc", p) for p in products)

                    if not reactants_have_boc and products_have_boc:
                        protection_found = True
                        print(f"Boc protection detected through FG analysis: {rsmi}")

                    # Check if Boc group is removed (deprotection)
                    if reactants_have_boc and not products_have_boc:
                        deprotection_found = True
                        print(f"Boc deprotection detected through FG analysis: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    # Consider a route to have Boc protection/deprotection strategy if Boc-protected molecules are found
    # This is more lenient than requiring both protection and deprotection reactions to be explicitly detected
    result = boc_protected_found
    print(f"Boc protection/deprotection strategy detected: {result}")
    print(
        f"Details: boc_protected_found={boc_protected_found}, protection_found={protection_found}, deprotection_found={deprotection_found}"
    )
    return result
