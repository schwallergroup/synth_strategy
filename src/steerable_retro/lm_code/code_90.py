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
    This function detects if the synthesis route preserves a cyclopropane ring
    from early starting materials through to the final product.
    """
    # Track if cyclopropane is preserved through the entire route
    preserved_through_route = [False]  # Using list for mutable reference in nested function

    def dfs_traverse(node, path=None):
        if path is None:
            path = []

        if node["type"] == "mol":
            # Check if this molecule contains a cyclopropane
            has_cyclopropane = checker.check_ring("cyclopropane", node["smiles"])

            # If this is the final product (depth 0), check if it has cyclopropane
            if len(path) == 0 and has_cyclopropane:
                print(f"Final product has cyclopropane: {node['smiles']}")
                path.append({"has_cyclopropane": has_cyclopropane, "smiles": node["smiles"]})
            elif has_cyclopropane:
                print(f"Intermediate at depth {len(path)} has cyclopropane: {node['smiles']}")
                path.append({"has_cyclopropane": has_cyclopropane, "smiles": node["smiles"]})
            else:
                print(f"Molecule at depth {len(path)} does NOT have cyclopropane: {node['smiles']}")
                path.append({"has_cyclopropane": False, "smiles": node["smiles"]})

            # If we've reached a starting material (leaf node)
            if node.get("in_stock", False) and len(node.get("children", [])) == 0:
                # Check if cyclopropane is preserved throughout this path
                cyclopropane_count = sum(1 for item in path if item["has_cyclopropane"])
                if cyclopropane_count == len(path) and len(path) > 1:
                    print(
                        f"Found complete path preserving cyclopropane from starting material to product!"
                    )
                    preserved_through_route[0] = True

                    # Print the path for debugging
                    for i, item in enumerate(path):
                        print(
                            f"Depth {i}: {'Has cyclopropane' if item['has_cyclopropane'] else 'No cyclopropane'} - {item['smiles']}"
                        )

        elif node["type"] == "reaction":
            # For reaction nodes, just continue traversal without modifying the path
            pass

        # Continue DFS traversal
        for child in node.get("children", []):
            # Create a copy of the path for this branch
            child_path = path.copy()
            dfs_traverse(child, child_path)

    # Start traversal from the root
    dfs_traverse(route)

    print(f"Cyclopropane preservation throughout route: {preserved_through_route[0]}")
    return preserved_through_route[0]
