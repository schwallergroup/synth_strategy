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
    This function detects a synthetic strategy involving thionation
    (conversion of amide to thioamide).
    """
    found_thionation = False

    def dfs_traverse(node):
        nonlocal found_thionation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                try:
                    rsmi = node["metadata"]["rsmi"]
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    print(f"Analyzing reaction: {rsmi}")

                    # Check if product contains thioamide
                    if checker.check_fg("Thioamide", product):
                        print(f"Product contains thioamide: {product}")

                        # Check if any reactant contains amide but not thioamide
                        for reactant in reactants:
                            # Check for primary, secondary, or tertiary amide
                            has_primary_amide = checker.check_fg("Primary amide", reactant)
                            has_secondary_amide = checker.check_fg("Secondary amide", reactant)
                            has_tertiary_amide = checker.check_fg("Tertiary amide", reactant)

                            has_amide = (
                                has_primary_amide or has_secondary_amide or has_tertiary_amide
                            )
                            has_thioamide = checker.check_fg("Thioamide", reactant)

                            if has_amide and not has_thioamide:
                                print(f"Reactant contains amide but not thioamide: {reactant}")
                                print(
                                    f"Primary amide: {has_primary_amide}, Secondary amide: {has_secondary_amide}, Tertiary amide: {has_tertiary_amide}"
                                )

                                # This is a thionation reaction
                                found_thionation = True
                                print("Found thionation of amide to thioamide")
                                break
                except Exception as e:
                    print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return found_thionation
