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
    This function detects if the synthesis route preserves both fluorine atoms
    and a cyano group throughout multiple steps.
    """
    # Track paths with preserved functional groups
    preservation_paths = []

    def dfs_traverse(node, current_path=None, depth=0):
        if current_path is None:
            current_path = []

        # For molecule nodes, check for fluorine and cyano groups
        if node["type"] == "mol" and node.get("smiles"):
            mol_smiles = node["smiles"]
            has_fluorine = checker.check_fg("Trifluoro group", mol_smiles) or checker.check_fg(
                "Aromatic halide", mol_smiles
            )
            has_cyano = checker.check_fg("Nitrile", mol_smiles)

            # Create or update the current path information
            current_info = {
                "depth": depth,
                "has_fluorine": has_fluorine,
                "has_cyano": has_cyano,
                "smiles": mol_smiles,
            }

            # Add to the current path
            new_path = current_path + [current_info]

            # If this is a leaf node (starting material), save the path
            if node.get("in_stock", False) or not node.get("children"):
                if len(new_path) > 1:  # Only save paths with multiple nodes
                    preservation_paths.append(new_path)

            # Continue traversal with updated path
            for child in node.get("children", []):
                dfs_traverse(child, new_path, depth + 1)

        # For reaction nodes, just continue traversal
        elif node["type"] == "reaction":
            for child in node.get("children", []):
                dfs_traverse(child, current_path, depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check each path for preservation of both groups for at least 3 consecutive steps
    for path in preservation_paths:
        # Sort by depth to ensure correct order (target to starting material)
        path.sort(key=lambda x: x["depth"])

        # Count consecutive steps with both groups preserved
        consecutive_count = 0
        max_consecutive = 0

        for node_info in path:
            if node_info["has_fluorine"] and node_info["has_cyano"]:
                consecutive_count += 1
                max_consecutive = max(max_consecutive, consecutive_count)
            else:
                consecutive_count = 0

        if max_consecutive >= 3:
            print(f"Found preservation path with {max_consecutive} consecutive steps")
            return True

    return False
