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
    Detects if a pyrimidine core is maintained throughout the synthesis.
    Only checks intermediates and the final product (non-in_stock molecules).
    Starting materials (in_stock=True) are excluded from this check.

    Args:
        route: A synthesis route dictionary following the SynthesisRoute schema

    Returns:
        bool: True if all non-starting material molecules contain a pyrimidine core
    """
    all_intermediates_have_pyrimidine = True

    def dfs_traverse(node):
        nonlocal all_intermediates_have_pyrimidine

        if node["type"] == "mol" and "smiles" in node:
            # Only check molecules that are not starting materials
            if not node.get("in_stock", False):
                try:
                    pyrimidine_present = checker.check_ring("pyrimidine", node["smiles"])
                    if not pyrimidine_present:
                        all_intermediates_have_pyrimidine = False
                        print(f"Intermediate without pyrimidine core: {node['smiles']}")
                    else:
                        print(f"Intermediate with pyrimidine core: {node['smiles']}")
                except Exception as e:
                    print(f"Error checking pyrimidine in {node['smiles']}: {e}")
                    all_intermediates_have_pyrimidine = False
            else:
                print(f"Skipping starting material: {node['smiles']}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Pyrimidine core maintained throughout intermediates: {all_intermediates_have_pyrimidine}"
    )
    return all_intermediates_have_pyrimidine
