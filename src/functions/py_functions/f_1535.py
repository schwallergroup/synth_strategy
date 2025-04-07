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
    Detects a synthetic strategy involving silyl protection of hydroxyl groups,
    particularly focusing on TBDMS (tert-butyldimethylsilyl) protection.
    """
    # Initialize tracking variables
    protection_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal protection_reactions

        if node["type"] == "reaction":
            try:
                # Extract reactants and products
                rsmi = node["metadata"]["rsmi"]

                # Check for silyl protection reaction
                if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                    protection_reactions += 1
                    print(
                        f"Detected silyl protection reaction at depth {depth}: {rsmi}"
                    )

                # Alternative check for TMS protection specifically
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if any reactant has alcohol and any product has silyl ether
                if product_smiles and any(reactants_smiles):
                    has_alcohol_reactant = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        or checker.check_fg("Aromatic alcohol", r)
                        for r in reactants_smiles
                    )

                    has_silyl_product = checker.check_fg(
                        "TMS ether protective group", product_smiles
                    ) or checker.check_fg("Silyl protective group", product_smiles)

                    if has_alcohol_reactant and has_silyl_product:
                        protection_reactions += 1
                        print(
                            f"Detected silyl protection (FG check) at depth {depth}: {rsmi}"
                        )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least one silyl protection reaction is detected
    strategy_detected = protection_reactions > 0
    print(
        f"Silyl protection strategy detected: {strategy_detected} (found {protection_reactions} reactions)"
    )
    return strategy_detected
