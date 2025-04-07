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
    This function detects if reductive amination is used to introduce amine functionality.
    It looks for aldehyde/ketone + amine -> amine transformation.
    """
    reductive_amination_found = False

    def dfs_traverse(node, depth=0):
        nonlocal reductive_amination_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]

            # Check if this is a reductive amination reaction directly
            if checker.check_reaction("reductive amination with aldehyde", rsmi):
                print(f"Found reductive amination with aldehyde at depth {depth}")
                reductive_amination_found = True
                return

            if checker.check_reaction("reductive amination with ketone", rsmi):
                print(f"Found reductive amination with ketone at depth {depth}")
                reductive_amination_found = True
                return

            if checker.check_reaction("reductive amination with alcohol", rsmi):
                print(f"Found reductive amination with alcohol at depth {depth}")
                reductive_amination_found = True
                return

            # If direct reaction check fails, check for the pattern manually
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            has_aldehyde = False
            has_ketone = False
            has_amine = False

            # Check reactants for required functional groups
            for reactant in reactants:
                if checker.check_fg("Aldehyde", reactant):
                    print(f"Found aldehyde in reactant: {reactant}")
                    has_aldehyde = True
                if checker.check_fg("Ketone", reactant):
                    print(f"Found ketone in reactant: {reactant}")
                    has_ketone = True
                if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                    "Secondary amine", reactant
                ):
                    print(f"Found amine in reactant: {reactant}")
                    has_amine = True

            # Check if product has amine but no carbonyl
            if (has_aldehyde or has_ketone) and has_amine:
                if (
                    checker.check_fg("Primary amine", product)
                    or checker.check_fg("Secondary amine", product)
                    or checker.check_fg("Tertiary amine", product)
                ):
                    # Check that product doesn't have the original carbonyl
                    if not (
                        checker.check_fg("Aldehyde", product) or checker.check_fg("Ketone", product)
                    ):
                        print(f"Found pattern consistent with reductive amination at depth {depth}")
                        print(f"Reactants: {reactants}")
                        print(f"Product: {product}")
                        reductive_amination_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return reductive_amination_found
