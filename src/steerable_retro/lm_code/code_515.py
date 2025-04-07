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
    Detects if the final step (depth 0) involves reduction of a nitro group to an amine.
    """
    print("Starting analysis to detect nitro reduction in final step")
    final_step_has_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_has_nitro_reduction

        print(f"Traversing node of type {node['type']} at depth {depth}")

        # Check if this is a reaction node
        if node["type"] == "reaction":
            print(f"Found reaction node at depth {depth}")

            # Check if this is the final step (depth 0)
            if depth == 0:
                print("Analyzing final step reaction (depth 0)")

                # Get reaction SMILES
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    print("No reaction SMILES found in metadata")
                else:
                    print(f"Reaction SMILES: {rsmi}")

                    # First check if this is a nitro reduction reaction directly
                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        print("Confirmed nitro reduction reaction via reaction checker")
                        final_step_has_nitro_reduction = True
                        return

                    # If specific reaction check fails, do more detailed analysis
                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        print(f"Reactants: {reactants}")
                        print(f"Product: {product}")

                        # Check if any reactant contains nitro group
                        has_nitro_reactant = any(
                            checker.check_fg("Nitro group", r) for r in reactants
                        )

                        # Check if product contains amine group
                        has_amine_product = (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                            or checker.check_fg("Aniline", product)
                        )

                        print(f"Has nitro in reactants: {has_nitro_reactant}")
                        print(f"Has amine in product: {has_amine_product}")

                        # If we have both nitro in reactants and amine in product, it's likely a reduction
                        if has_nitro_reactant and has_amine_product:
                            # Additional check: make sure the nitro group is actually gone in the product
                            if not checker.check_fg("Nitro group", product):
                                print("Found nitro reduction to amine in final step")
                                final_step_has_nitro_reduction = True
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # If we didn't find a reaction at depth 0, try to find the reaction closest to the target molecule
    if not final_step_has_nitro_reduction and route["type"] == "mol":
        print("No reaction found at depth 0, checking first child reaction")
        for child in route.get("children", []):
            if child["type"] == "reaction":
                print("Found first child reaction, analyzing it as final step")
                rsmi = child.get("metadata", {}).get("rsmi", "")
                if rsmi:
                    print(f"Reaction SMILES: {rsmi}")

                    # Check if this is a nitro reduction
                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        print("Confirmed nitro reduction reaction via reaction checker")
                        final_step_has_nitro_reduction = True
                        break

                    try:
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check if any reactant contains nitro group
                        has_nitro_reactant = any(
                            checker.check_fg("Nitro group", r) for r in reactants
                        )

                        # Check if product contains amine group
                        has_amine_product = (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                            or checker.check_fg("Aniline", product)
                        )

                        # If we have both nitro in reactants and amine in product, it's likely a reduction
                        if has_nitro_reactant and has_amine_product:
                            # Additional check: make sure the nitro group is actually gone in the product
                            if not checker.check_fg("Nitro group", product):
                                print("Found nitro reduction to amine in first reaction")
                                final_step_has_nitro_reduction = True
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")
                break

    print(f"Final result: {final_step_has_nitro_reduction}")
    return final_step_has_nitro_reduction
