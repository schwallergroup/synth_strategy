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
    Detects if the synthesis route contains reductive amination in early steps (high depth).
    """
    has_early_reductive_amination = False

    def dfs_traverse(node, depth=0):
        nonlocal has_early_reductive_amination

        if node["type"] == "reaction" and depth >= 2:  # Early steps have depth >= 2
            try:
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Direct check for reductive amination reaction
                if checker.check_reaction("reductive amination with aldehyde", rsmi):
                    print(f"Detected reductive amination with aldehyde at depth {depth}")
                    has_early_reductive_amination = True
                elif checker.check_reaction("reductive amination with ketone", rsmi):
                    print(f"Detected reductive amination with ketone at depth {depth}")
                    has_early_reductive_amination = True
                elif checker.check_reaction("reductive amination with alcohol", rsmi):
                    print(f"Detected reductive amination with alcohol at depth {depth}")
                    has_early_reductive_amination = True
                else:
                    # Fallback check for reductive amination pattern
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for carbonyl compounds and amines in reactants
                    has_carbonyl = any(
                        checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                        for r in reactants
                    )

                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        for r in reactants
                    )

                    # Check for secondary or tertiary amine in product
                    has_product_amine = checker.check_fg(
                        "Secondary amine", product
                    ) or checker.check_fg("Tertiary amine", product)

                    if has_carbonyl and has_amine and has_product_amine:
                        print(f"Detected potential reductive amination pattern at depth {depth}")
                        has_early_reductive_amination = True
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: {has_early_reductive_amination}")
    return has_early_reductive_amination
