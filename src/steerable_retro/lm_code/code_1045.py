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
    This function detects a late-stage alcohol protection with a silyl group.
    """
    # Track if we found the silyl protection reaction
    found_silyl_protection = False

    def manual_silyl_protection_check(reactants, product):
        """
        Manually check if the reaction is a silyl protection by examining
        reactants and products for characteristic functional groups.
        """
        # Check if any reactant has a silyl group (likely a silyl chloride or similar)
        has_silyl_reactant = any(
            checker.check_fg("Silyl protective group", r)
            or "Si" in r  # Simple check for silicon atom
            for r in reactants
        )

        # Check if any reactant has an alcohol
        alcohol_types = [
            "Primary alcohol",
            "Secondary alcohol",
            "Tertiary alcohol",
            "Aromatic alcohol",
        ]
        has_alcohol_reactant = any(
            any(checker.check_fg(alcohol_type, r) for alcohol_type in alcohol_types)
            for r in reactants
        )

        # Check if product has a silyl ether
        has_silyl_ether_product = (
            checker.check_fg("TMS ether protective group", product)
            or checker.check_fg("Silyl protective group", product)
            or (
                # Check for Si-O pattern in product but not in reactants
                "Si" in product
                and "O" in product
                and not any(("Si" in r and "O" in r) for r in reactants)
            )
        )

        # For a silyl protection, we need a silyl reagent, an alcohol, and the product should have a silyl ether
        return has_silyl_reactant and has_alcohol_reactant and has_silyl_ether_product

    def dfs_traverse(node, depth=0):
        nonlocal found_silyl_protection

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}, RSMI: {rsmi}")

            # Check if this is a late-stage reaction (depth 0-3)
            if depth <= 3:
                print(f"This is a late-stage reaction (depth {depth})")

                # Check if this is a silyl protection reaction
                is_silyl_protection = checker.check_reaction(
                    "Alcohol protection with silyl ethers", rsmi
                )
                print(f"Is silyl protection reaction (from checker): {is_silyl_protection}")

                # If not detected by the checker, try manual detection
                if not is_silyl_protection:
                    is_silyl_protection = manual_silyl_protection_check(reactants, product)
                    print(
                        f"Is silyl protection reaction (from manual check): {is_silyl_protection}"
                    )

                if is_silyl_protection:
                    # Verify reactant has alcohol and product has silyl ether
                    alcohol_types = [
                        "Primary alcohol",
                        "Secondary alcohol",
                        "Tertiary alcohol",
                        "Aromatic alcohol",
                    ]
                    has_alcohol_reactant = any(
                        any(checker.check_fg(alcohol_type, r) for alcohol_type in alcohol_types)
                        for r in reactants
                    )

                    has_silyl_ether_product = (
                        checker.check_fg("TMS ether protective group", product)
                        or checker.check_fg("Silyl protective group", product)
                        or (
                            # Check for Si-O pattern in product but not in reactants
                            "Si" in product
                            and "O" in product
                            and not any(("Si" in r and "O" in r) for r in reactants)
                        )
                    )

                    print(f"Has alcohol reactant: {has_alcohol_reactant}")
                    print(f"Has silyl ether product: {has_silyl_ether_product}")

                    if has_alcohol_reactant and has_silyl_ether_product:
                        found_silyl_protection = True
                        print(f"Found late-stage silyl protection reaction at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_silyl_protection
