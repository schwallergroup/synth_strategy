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
    Detects if heterocycles (pyrazole, pyridine, oxetane) are preserved throughout the synthesis.
    """
    # List of heterocycles to track
    heterocycles = [
        "pyrazole",
        "pyridine",
        "oxetane",
        "furan",
        "thiophene",
        "imidazole",
        "thiazole",
        "oxazole",
        "isoxazole",
        "pyrimidine",
        "piperidine",
        "morpholine",
    ]

    # Track heterocycles at each depth and their atom indices
    heterocycles_at_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                # Check for each heterocycle
                for heterocycle in heterocycles:
                    # Check if molecule contains the heterocycle
                    has_heterocycle = checker.check_ring(heterocycle, smiles)

                    # If heterocycle is found, get its atom indices
                    if has_heterocycle:
                        # Get atom indices for the heterocycle
                        indices = checker.get_ring_atom_indices(heterocycle, smiles)

                        # Initialize depth entry if not exists
                        if depth not in heterocycles_at_depth:
                            heterocycles_at_depth[depth] = {}

                        # Initialize heterocycle entry if not exists
                        if heterocycle not in heterocycles_at_depth[depth]:
                            heterocycles_at_depth[depth][heterocycle] = 0

                        # Increment count of this heterocycle at this depth
                        heterocycles_at_depth[depth][heterocycle] += len(indices)

                        print(
                            f"Depth {depth}: Found {heterocycle}, count: {heterocycles_at_depth[depth][heterocycle]}"
                        )

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if heterocycles are preserved throughout
    if not heterocycles_at_depth:
        print("No heterocycles found in the route")
        return False

    max_depth = max(heterocycles_at_depth.keys())
    print(f"Max depth: {max_depth}")
    print(f"Heterocycles at each depth: {heterocycles_at_depth}")

    # Get all depths where molecules were found
    depths = sorted(heterocycles_at_depth.keys())

    # Check if at least two heterocycles are preserved from beginning to end
    preserved_count = 0
    for heterocycle in heterocycles:
        preserved = True
        # Check if this heterocycle exists at depth 0 (final product)
        if (
            0 not in heterocycles_at_depth
            or heterocycle not in heterocycles_at_depth[0]
            or heterocycles_at_depth[0][heterocycle] == 0
        ):
            preserved = False
            continue

        # Check if this heterocycle exists at all depths where molecules were found
        for depth in depths:
            if (
                heterocycle not in heterocycles_at_depth[depth]
                or heterocycles_at_depth[depth][heterocycle] == 0
            ):
                preserved = False
                break

        if preserved:
            print(f"Heterocycle {heterocycle} is preserved throughout the synthesis")
            preserved_count += 1

            # Early return if we've found at least two preserved heterocycles
            if preserved_count >= 2:
                print(f"Found {preserved_count} preserved heterocycles, returning True")
                return True

    print(f"Preserved heterocycles count: {preserved_count}")
    return preserved_count >= 2
