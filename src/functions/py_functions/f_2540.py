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
    Detects a synthetic strategy involving Boc protection of an amine
    that is maintained throughout the synthesis.
    """
    boc_protection_depth = None
    boc_deprotection_depth = None
    max_depth = 0
    boc_maintained = True
    boc_actually_added = False
    final_product_has_boc = False
    starting_material_has_boc = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_depth, boc_deprotection_depth, max_depth
        nonlocal boc_maintained, boc_actually_added, final_product_has_boc, starting_material_has_boc

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[2]

            # Check for Boc protection reaction
            if (
                checker.check_reaction("Boc amine protection", rsmi)
                or checker.check_reaction("Boc amine protection explicit", rsmi)
                or checker.check_reaction(
                    "Boc amine protection with Boc anhydride", rsmi
                )
                or checker.check_reaction("Boc amine protection of primary amine", rsmi)
                or checker.check_reaction(
                    "Boc amine protection of secondary amine", rsmi
                )
                or checker.check_reaction("Boc amine protection (ethyl Boc)", rsmi)
            ):

                # Verify Boc was actually added (not in reactants, but in product)
                reactants_have_boc = any(
                    checker.check_fg("Boc", reactant) for reactant in reactants
                )
                product_has_boc = checker.check_fg("Boc", product)

                if product_has_boc and not reactants_have_boc:
                    boc_protection_depth = depth
                    boc_actually_added = True
                    print(f"Detected Boc protection at depth {depth} - Boc was added")

            # Check if this is a Boc deprotection reaction
            if (
                checker.check_reaction("Boc amine deprotection", rsmi)
                or checker.check_reaction("Boc amine deprotection of guanidine", rsmi)
                or checker.check_reaction("Boc amine deprotection to NH-NH2", rsmi)
                or checker.check_reaction("Tert-butyl deprotection of amine", rsmi)
            ):

                # Verify Boc was actually removed
                reactants_have_boc = any(
                    checker.check_fg("Boc", reactant) for reactant in reactants
                )
                product_has_boc = checker.check_fg("Boc", product)

                if reactants_have_boc and not product_has_boc:
                    boc_deprotection_depth = depth
                    print(f"Detected Boc deprotection at depth {depth}")

            # Check if Boc group is maintained in subsequent reactions
            if (
                boc_protection_depth is not None
                and boc_deprotection_depth is None
                and depth > boc_protection_depth
            ):
                # Check if product still has Boc group
                if not checker.check_fg("Boc", product):
                    # Verify this isn't a deprotection reaction
                    if not any(
                        checker.check_reaction(rxn, rsmi)
                        for rxn in [
                            "Boc amine deprotection",
                            "Boc amine deprotection of guanidine",
                            "Boc amine deprotection to NH-NH2",
                            "Tert-butyl deprotection of amine",
                        ]
                    ):
                        boc_maintained = False
                        print(f"Boc group unexpectedly removed at depth {depth}")

        elif node["type"] == "mol":
            if depth == 0:  # This is the final product
                final_product_has_boc = checker.check_fg("Boc", node["smiles"])
                print(f"Final product contains a Boc group: {final_product_has_boc}")

            # Check if this is a starting material with Boc
            if node.get("in_stock", False) and checker.check_fg("Boc", node["smiles"]):
                starting_material_has_boc = True
                print(f"Starting material at depth {depth} contains a Boc group")

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # A proper Boc protection strategy can be:
    # 1. Add a Boc group during synthesis and maintain it
    # 2. Start with a Boc-protected compound and maintain it throughout

    # Case 1: Boc added during synthesis
    strategy_case1 = (
        boc_protection_depth is not None and boc_actually_added and boc_maintained
    )

    # Case 2: Boc present in starting material and maintained
    strategy_case2 = (
        starting_material_has_boc and final_product_has_boc and boc_maintained
    )

    strategy_present = strategy_case1 or strategy_case2

    if strategy_present:
        if strategy_case1:
            print(
                f"Detected Boc protection strategy: protected at depth {boc_protection_depth} and maintained throughout synthesis"
            )
        else:
            print(
                "Detected Boc protection strategy: Boc present in starting material and maintained throughout synthesis"
            )

        if boc_deprotection_depth is not None:
            print(f"Boc group was deprotected at depth {boc_deprotection_depth}")
        else:
            print("Boc group was maintained in the final product")

    return strategy_present
