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
    Detects if a morpholine ring is preserved throughout the synthesis.
    """
    steps_with_morpholine = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal steps_with_morpholine, total_steps

        if node["type"] == "mol" and node.get("smiles"):
            try:
                # Only count non-stock molecules (intermediates and target)
                if not node.get("in_stock", False):
                    # Use the checker function to detect morpholine
                    if checker.check_ring("morpholine", node["smiles"]):
                        steps_with_morpholine += 1
                        print(f"Found morpholine in: {node['smiles']}")
                    total_steps += 1
            except Exception as e:
                print(
                    f"Error processing SMILES in morpholine detection: {node['smiles']}, Error: {e}"
                )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Steps with morpholine: {steps_with_morpholine}, Total steps: {total_steps}")
    # Return True if morpholine is present in most molecules (>75%)
    return steps_with_morpholine > 0 and steps_with_morpholine / max(1, total_steps) > 0.75
