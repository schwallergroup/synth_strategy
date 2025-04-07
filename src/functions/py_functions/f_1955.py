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
    This function detects an ester to alcohol reduction in the synthetic route.
    """
    found_ester_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal found_ester_reduction

        # Check for reaction nodes
        if node["type"] == "reaction":
            # Get reaction SMILES
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # First check if this is a known ester reduction reaction
                if checker.check_reaction(
                    "Reduction of ester to primary alcohol", rsmi
                ):
                    print(
                        f"Found ester to alcohol reduction via reaction check: {rsmi}"
                    )
                    found_ester_reduction = True
                    return

                # If reaction check fails, check for functional group transformation
                # Check if any reactant contains an ester group
                has_ester_reactant = any(
                    checker.check_fg("Ester", r) for r in reactants
                )
                # Check if product contains a primary alcohol
                has_primary_alcohol_product = checker.check_fg(
                    "Primary alcohol", product
                )

                if has_ester_reactant and has_primary_alcohol_product:
                    # Additional check to ensure this is a reduction reaction
                    # In a reduction, we should not see primary alcohol in reactants
                    reactants_with_primary_alcohol = sum(
                        1 for r in reactants if checker.check_fg("Primary alcohol", r)
                    )

                    # If product has more primary alcohols than reactants, it's likely a reduction
                    if reactants_with_primary_alcohol < has_primary_alcohol_product:
                        print(
                            f"Found ester to alcohol reduction via functional group check:"
                        )
                        print(f"  Reactants: {reactants}")
                        print(f"  Product: {product}")
                        found_ester_reduction = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_ester_reduction
