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
    Detects late-stage sulfonamide formation in the synthetic route.
    """
    sulfonamide_found = False
    late_stage = False
    max_depth_found = 0
    sulfonamide_depth = float("inf")

    def dfs(node, depth=0):
        nonlocal sulfonamide_found, late_stage, max_depth_found, sulfonamide_depth

        # Update max depth to determine late vs early stage
        max_depth_found = max(max_depth_found, depth)

        if node["type"] == "mol":
            # Check if molecule contains sulfonamide group
            if checker.check_fg("Sulfonamide", node["smiles"]):
                print(f"Found sulfonamide group in molecule at depth {depth}: {node['smiles']}")
                sulfonamide_found = True
                sulfonamide_depth = min(sulfonamide_depth, depth)

        elif node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]

                # Check for sulfonamide formation reactions
                if checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                ):
                    print(f"Found sulfonamide formation reaction at depth {depth}: {rsmi}")
                    sulfonamide_found = True
                    sulfonamide_depth = min(sulfonamide_depth, depth)
            except Exception as e:
                print(f"Error in sulfonamide check: {e}")

        # Continue DFS traversal
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)

    # Determine if sulfonamide formation is late-stage
    # Late stage is typically in the first third of the synthesis depth
    if max_depth_found > 0 and sulfonamide_found:
        if sulfonamide_depth <= 2:  # Absolute measure: depth <= 2 is late stage
            late_stage = True
        elif max_depth_found > 6:  # Relative measure for deeper trees
            late_stage = sulfonamide_depth <= max_depth_found // 3
        else:
            late_stage = sulfonamide_depth <= 2  # Default for medium-sized trees

    print(
        f"Late-stage sulfonamide formation detected: {late_stage} (sulfonamide found: {sulfonamide_found}, depth: {sulfonamide_depth}, max depth: {max_depth_found})"
    )
    return late_stage
