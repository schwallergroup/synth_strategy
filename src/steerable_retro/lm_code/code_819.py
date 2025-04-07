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
    This function detects if the route contains formation of a tertiary alcohol from a methylene group.
    """
    found_tertiary_alcohol = False

    def dfs_traverse(node):
        nonlocal found_tertiary_alcohol

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if tertiary alcohol is in product but not in reactants
            if checker.check_fg("Tertiary alcohol", product_part) and not checker.check_fg(
                "Tertiary alcohol", reactants_part
            ):
                # Check for relevant reaction types that form tertiary alcohols
                if (
                    checker.check_reaction("Grignard from ketone to alcohol", rsmi)
                    or checker.check_reaction("Grignard from aldehyde to alcohol", rsmi)
                    or checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                ):

                    # Check for methylene or vinyl group in reactants
                    if checker.check_fg("Vinyl", reactants_part) or checker.check_fg(
                        "Ethylene", reactants_part
                    ):
                        print(f"Found tertiary alcohol formation reaction: {rsmi}")
                        found_tertiary_alcohol = True

                # Additional check for other reaction types that might form tertiary alcohols
                elif checker.check_fg("Alkyne", reactants_part) and "OH" in product_part:
                    print(f"Found tertiary alcohol formation from alkyne: {rsmi}")
                    found_tertiary_alcohol = True

                # Check for addition to carbonyl compounds
                elif (
                    checker.check_fg("Aldehyde", reactants_part)
                    or checker.check_fg("Ketone", reactants_part)
                ) and "OH" in product_part:
                    print(f"Found tertiary alcohol formation from carbonyl: {rsmi}")
                    found_tertiary_alcohol = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_tertiary_alcohol
