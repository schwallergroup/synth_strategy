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
    Detects the use of Boc protection of amines in the synthesis.
    """
    boc_protection_found = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_found

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            try:
                rsmi = node["metadata"]["rsmi"]

                # Check for both Boc protection and deprotection reactions
                boc_reactions = [
                    # Protection reactions
                    "Boc amine protection",
                    "Boc amine protection explicit",
                    "Boc amine protection with Boc anhydride",
                    "Boc amine protection (ethyl Boc)",
                    "Boc amine protection of secondary amine",
                    "Boc amine protection of primary amine",
                    # Deprotection reactions
                    "Boc amine deprotection",
                    "Boc amine deprotection of guanidine",
                    "Boc amine deprotection to NH-NH2",
                    "Tert-butyl deprotection of amine",
                ]

                for reaction_type in boc_reactions:
                    if checker.check_reaction(reaction_type, rsmi):
                        print(f"Found Boc reaction at depth {depth}: {reaction_type}")
                        boc_protection_found = True
                        # Continue traversal to find all instances

                # If no direct reaction match, check for Boc group addition/removal
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for Boc addition (protection)
                product_has_boc = checker.check_fg("Carbamic ester", product)
                reactants_have_boc = any(
                    checker.check_fg("Carbamic ester", r) for r in reactants
                )

                if product_has_boc and not reactants_have_boc:
                    print(
                        f"Found Boc protection via functional group analysis at depth {depth}"
                    )
                    boc_protection_found = True

                # Check for Boc removal (deprotection)
                product_has_amine = checker.check_fg(
                    "Primary amine", product
                ) or checker.check_fg("Secondary amine", product)
                if product_has_amine and reactants_have_boc:
                    print(
                        f"Found Boc deprotection via functional group analysis at depth {depth}"
                    )
                    boc_protection_found = True

            except Exception as e:
                print(f"Error analyzing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    if boc_protection_found:
        print("Detected Boc protection strategy")
    else:
        print("No Boc protection strategy detected")

    return boc_protection_found
