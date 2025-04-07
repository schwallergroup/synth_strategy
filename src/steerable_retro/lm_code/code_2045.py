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
    Detects synthesis routes involving oxidation of secondary alcohols to ketones.
    In retrosynthetic analysis, this means detecting reduction of ketones to secondary alcohols.
    """
    oxidation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal oxidation_detected

        if node["type"] == "reaction":
            try:
                rsmi = node.get("metadata", {}).get("rsmi", "")
                if not rsmi:
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Depth {depth}: Analyzing reaction: {rsmi}")

                # Check for known reduction reactions (oxidation in forward direction)
                reduction_reactions = [
                    "Reduction of ketone to secondary alcohol",
                    "Reduction of aldehydes and ketones to alcohols",
                ]

                for rxn_type in reduction_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Depth {depth}: Detected {rxn_type} reaction")
                        oxidation_detected = True
                        return

                # Check for oxidation reactions (which we interpret in reverse for retrosynthesis)
                oxidation_reactions = [
                    "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones",
                    "Oxidation of alcohol to carboxylic acid",
                ]

                for rxn_type in oxidation_reactions:
                    if checker.check_reaction(rxn_type, rsmi):
                        print(f"Depth {depth}: Detected {rxn_type} reaction")

                        # In retrosynthesis, ketone should be in reactants, alcohol in product
                        has_ketone = any(
                            checker.check_fg("Ketone", reactant) for reactant in reactants_smiles
                        )
                        has_sec_alcohol = checker.check_fg("Secondary alcohol", product_smiles)

                        if has_ketone and has_sec_alcohol:
                            print(
                                f"Depth {depth}: Confirmed ketone reduction to secondary alcohol (oxidation in forward)"
                            )
                            oxidation_detected = True
                            return

                # Alternative check for cases where the reaction might not be directly recognized
                # Check if reactants contain ketone and product contains secondary alcohol
                has_ketone_reactants = any(
                    checker.check_fg("Ketone", reactant) for reactant in reactants_smiles
                )
                has_sec_alcohol_product = checker.check_fg("Secondary alcohol", product_smiles)

                # In retrosynthesis: ketone in reactants, alcohol in product
                if has_ketone_reactants and has_sec_alcohol_product:
                    print(
                        f"Depth {depth}: Detected potential ketone reduction to secondary alcohol"
                    )
                    oxidation_detected = True
                    return

                # Special case for reactions like at Depth 7
                # Look for carbonyl to alcohol conversion using atom mapping
                if "[C" in rsmi and "=[O" in rsmi.split(">")[0] and "[OH" in rsmi.split(">")[-1]:
                    # Extract atom mappings to check if the same carbon is involved
                    for reactant in reactants_smiles:
                        if "=[O:" in reactant:
                            # Extract the atom mapping number for the oxygen
                            o_mapping = reactant.split("=[O:")[1].split("]")[0]
                            # Check if the same mapping appears as hydroxyl in product
                            if (
                                f"[OH:{o_mapping}]" in product_smiles
                                or f"OH:{o_mapping}" in product_smiles
                            ):
                                print(
                                    f"Depth {depth}: Detected carbonyl reduction to alcohol based on atom mapping"
                                )
                                oxidation_detected = True
                                return

            except Exception as e:
                print(f"Error in alcohol oxidation detection: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Final result: oxidation_detected = {oxidation_detected}")

    return oxidation_detected
