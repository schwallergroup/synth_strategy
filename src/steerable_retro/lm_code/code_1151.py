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
    Detects a synthetic strategy where a trifluoromethyl group is preserved
    throughout the synthesis.
    """
    # Track if CF3 is present in all molecules
    cf3_present_in_all = True

    def dfs_traverse(node, depth=0):
        nonlocal cf3_present_in_all

        if node["type"] == "mol" and not node.get("in_stock", False):
            # Check for CF3 group in non-starting materials
            print(f"Checking molecule at depth {depth}: {node['smiles']}")

            # Use the checker function to detect trifluoromethyl group
            has_cf3 = checker.check_fg("Trifluoro group", node["smiles"])

            if not has_cf3:
                cf3_present_in_all = False
                print(f"Molecule without CF3 group found: {node['smiles']}")
            else:
                print(f"CF3 group found in molecule: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"CF3 preservation strategy detected: {cf3_present_in_all}")
    return cf3_present_in_all
