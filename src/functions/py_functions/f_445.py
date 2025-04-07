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
    Detects if the synthesis route contains and preserves a trifluoromethyl group.
    """
    # Check if the final product contains a trifluoromethyl group
    if route["type"] != "mol" or not checker.check_fg(
        "Trifluoro group", route["smiles"]
    ):
        return False

    # Track if we find a trifluoromethyl group in any starting material
    trifluoro_in_starting_material = False

    def dfs(node, depth=0):
        nonlocal trifluoro_in_starting_material

        if node["type"] == "mol" and node.get("in_stock", False) and node["smiles"]:
            if checker.check_fg("Trifluoro group", node["smiles"]):
                trifluoro_in_starting_material = True
                print(f"Found trifluoromethyl in starting material: {node['smiles']}")

        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)

    return trifluoro_in_starting_material
