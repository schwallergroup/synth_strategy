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
    Detects a synthesis strategy where certain functional groups (like methoxy)
    are preserved throughout the synthesis.
    """
    # Track molecules with methoxy groups at each depth
    molecules_with_methoxy = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            print(f"Checking molecule at depth {depth}: {mol_smiles}")

            # Check for methoxy group using the checker
            has_methoxy = checker.check_fg("Ether", mol_smiles)

            if has_methoxy:
                print(f"Found methoxy group at depth {depth}")
                # Store molecule with methoxy at this depth
                if depth not in molecules_with_methoxy:
                    molecules_with_methoxy[depth] = []
                molecules_with_methoxy[depth].append(mol_smiles)

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Count depths where methoxy was found
    depths_with_methoxy = list(molecules_with_methoxy.keys())
    depths_with_methoxy.sort()  # Sort depths for analysis

    print(f"Found methoxy groups at depths: {depths_with_methoxy}")

    # Check if methoxy is preserved across at least 3 steps (depths)
    # We need to find at least 3 different depths with methoxy groups
    result = len(depths_with_methoxy) >= 3

    # Additional check: ensure the depths are reasonably distributed
    # This helps confirm the methoxy group is preserved throughout the synthesis
    if result and len(depths_with_methoxy) >= 3:
        # Check if the methoxy appears in early, middle, and late stages
        # by dividing the depth range into three segments
        max_depth = max(depths_with_methoxy)
        if max_depth >= 2:  # Ensure we have enough depth range
            early_stage = set(range(max_depth * 2 // 3, max_depth + 1))
            middle_stage = set(range(max_depth // 3, max_depth * 2 // 3))
            late_stage = set(range(0, max_depth // 3))

            has_early = any(d in early_stage for d in depths_with_methoxy)
            has_middle = any(d in middle_stage for d in depths_with_methoxy)
            has_late = any(d in late_stage for d in depths_with_methoxy)

            # Require methoxy to appear in at least early and late stages
            result = has_early and has_late
            print(
                f"Distribution check: early={has_early}, middle={has_middle}, late={has_late}"
            )

    print(
        f"Preserved methoxy group strategy: {result}, found at {len(depths_with_methoxy)} depths"
    )
    return result
