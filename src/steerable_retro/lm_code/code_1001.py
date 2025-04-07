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
    Detects if the synthesis is centered around a heterocyclic structure
    that persists throughout the synthesis.
    """
    # Common heterocycles to check
    heterocycles = [
        "pyrimidine",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "triazole",
        "tetrazole",
        "furan",
        "thiophene",
        "pyrrole",
        "indole",
        "benzimidazole",
        "benzoxazole",
        "benzothiazole",
    ]

    # Track depths and molecules where heterocycles are found
    heterocycle_info = {}  # {heterocycle_type: {depth: smiles}}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]
            try:
                # Check for each heterocycle type
                for heterocycle in heterocycles:
                    if checker.check_ring(heterocycle, mol_smiles):
                        if heterocycle not in heterocycle_info:
                            heterocycle_info[heterocycle] = {}
                        heterocycle_info[heterocycle][depth] = mol_smiles
                        print(f"{heterocycle} detected at depth {depth} in molecule: {mol_smiles}")
            except Exception as e:
                print(f"Error checking heterocycle: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if any heterocycle persists through multiple depths
    for heterocycle, depths_dict in heterocycle_info.items():
        if len(depths_dict) >= 2:
            print(
                f"Heterocycle-based synthesis strategy detected: {heterocycle} persists through {len(depths_dict)} depths"
            )
            return True

    print("No heterocycle-based synthesis strategy detected")
    return False
