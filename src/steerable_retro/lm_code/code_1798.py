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
    This function detects a synthetic strategy involving reduction of a ketone to a secondary alcohol.
    """
    # Track if we found the ketone reduction
    found_ketone_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ketone_reduction

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reaction information
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a ketone reduction reaction
            is_ketone_reduction = checker.check_reaction(
                "Reduction of ketone to secondary alcohol", rsmi
            )

            if is_ketone_reduction:
                # In retrosynthesis, the product contains the ketone and reactants contain the alcohol
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                # Check for secondary alcohol in reactants and ketone in product
                has_secondary_alcohol = checker.check_fg("Secondary alcohol", reactants_part)
                has_ketone = checker.check_fg("Ketone", product_part)

                if has_ketone and has_secondary_alcohol:
                    print(f"Found ketone reduction at depth {depth}, reaction: {rsmi}")
                    found_ketone_reduction = True
                else:
                    print(
                        f"Reaction classified as ketone reduction but FG check failed at depth {depth}"
                    )
                    print(f"  Has ketone in product: {has_ketone}")
                    print(f"  Has secondary alcohol in reactant: {has_secondary_alcohol}")
            else:
                # Even if not classified as ketone reduction, check if it might be one
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]

                has_secondary_alcohol = checker.check_fg("Secondary alcohol", reactants_part)
                has_ketone = checker.check_fg("Ketone", product_part)

                if has_ketone and has_secondary_alcohol:
                    print(
                        f"Potential unlabeled ketone reduction at depth {depth}, reaction: {rsmi}"
                    )
                    found_ketone_reduction = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return found_ketone_reduction
