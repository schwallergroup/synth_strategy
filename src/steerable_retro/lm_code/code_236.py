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
    Detects if the synthesis route contains a Grignard reaction (ketone to tertiary alcohol).
    """
    grignard_found = False

    def dfs_traverse(node, depth=0):
        nonlocal grignard_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check specifically for Grignard reaction from ketone to alcohol
            if checker.check_reaction("Grignard from ketone to alcohol", rsmi):
                print(f"Grignard reaction (ketone to alcohol) found at depth {depth}")
                grignard_found = True

            # Alternative check for Grignard reactions
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains magnesium and another contains a ketone
            has_mg = any(checker.check_fg("Magnesium halide", r) for r in reactants)
            has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)

            # Check if product contains a tertiary alcohol
            if has_mg and has_ketone and product:
                if checker.check_fg("Tertiary alcohol", product):
                    print(f"Grignard pattern detected manually at depth {depth}")
                    grignard_found = True

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return grignard_found
