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
    Detects tosylation activation strategy in the synthetic route.
    """
    tosylation_found = False

    def dfs(node, depth=0):
        nonlocal tosylation_found

        if tosylation_found:
            return

        if node["type"] == "mol":
            # Check if molecule contains tosylate group
            if checker.check_fg("Tosylate", node["smiles"]):
                print(f"Found tosylate group in molecule: {node['smiles']}")
                tosylation_found = True

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]

                # Check for tosylation reaction
                if checker.check_reaction("Formation of Sulfonic Esters", rsmi):
                    print(f"Found tosylation reaction: {rsmi}")
                    tosylation_found = True
            except Exception as e:
                print(f"Error in tosylation check: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    print(f"Tosylation activation strategy detected: {tosylation_found}")
    return tosylation_found
