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
    Detects if the synthesis involves reduction of a ketone or aldehyde to a secondary alcohol.
    """
    carbonyl_reduction_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal carbonyl_reduction_detected

        if carbonyl_reduction_detected:
            return  # Stop traversal if already detected

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # First try to directly check for the specific reaction types
            if checker.check_reaction(
                "Reduction of ketone to secondary alcohol", rsmi
            ) or checker.check_reaction(
                "Reduction of aldehydes and ketones to alcohols", rsmi
            ):
                print(f"Detected carbonyl reduction reaction: {rsmi}")
                carbonyl_reduction_detected = True
                return

            # Alternative approach: manually check for carbonyl reduction
            try:
                reactants_part = rsmi.split(">")[0]
                product_part = rsmi.split(">")[-1]
                reactants = reactants_part.split(".")

                # Check if any reactant has a ketone or aldehyde and product has a secondary alcohol
                carbonyl_in_reactants = any(
                    checker.check_fg("Ketone", r)
                    or checker.check_fg("Aldehyde", r)
                    or checker.check_fg("Formaldehyde", r)
                    for r in reactants
                )
                sec_alcohol_in_product = checker.check_fg(
                    "Secondary alcohol", product_part
                )

                if carbonyl_in_reactants and sec_alcohol_in_product:
                    print(
                        f"Detected carbonyl group in reactants and secondary alcohol in product: {rsmi}"
                    )
                    carbonyl_reduction_detected = True

                # Special case: check for specific atom transformations
                # This handles cases where the functional group detection might miss some patterns
                if not carbonyl_reduction_detected:
                    for reactant in reactants:
                        if "[CH](=[O])" in reactant or "C(=[O])" in reactant:
                            if (
                                "[CH]([OH])" in product_part
                                or "C([OH])" in product_part
                            ):
                                print(
                                    f"Detected carbonyl to alcohol transformation based on atom patterns: {rsmi}"
                                )
                                carbonyl_reduction_detected = True
                                break

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return carbonyl_reduction_detected
