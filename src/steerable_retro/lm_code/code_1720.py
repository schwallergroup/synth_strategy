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
    Detects if the synthesis follows a linear strategy with sequential
    transformations on heterocyclic scaffolds.
    """
    # List of heterocyclic rings to check
    heterocycle_types = [
        "furan",
        "pyran",
        "pyrrole",
        "pyridine",
        "pyrazole",
        "imidazole",
        "oxazole",
        "thiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
        "triazole",
        "tetrazole",
        "indole",
        "quinoline",
        "isoquinoline",
        "purine",
        "thiophene",
        "benzoxazole",
        "benzothiazole",
        "benzimidazole",
        "isoxazole",
        "isothiazole",
        "oxadiazole",
        "thiadiazole",
    ]

    # Track the paths of heterocycle transformations
    paths = []
    current_path = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            # Check if this molecule contains a heterocycle
            has_heterocycle = False
            heterocycle_names = []

            for het_type in heterocycle_types:
                if checker.check_ring(het_type, mol_smiles):
                    has_heterocycle = True
                    heterocycle_names.append(het_type)

            if has_heterocycle:
                # Add this molecule to the current path
                current_path.append(
                    {"depth": depth, "smiles": mol_smiles, "heterocycles": heterocycle_names}
                )
                print(
                    f"Depth {depth}: Found molecule with heterocycles: {', '.join(heterocycle_names)}"
                )

            # If this is a leaf node (starting material) and we have a path with heterocycles
            if node.get("in_stock", False) and len(current_path) > 0:
                # Save this path
                paths.append(list(current_path))
                print(f"Completed path with {len(current_path)} heterocycle-containing molecules")

        # Process children (in retrosynthetic direction)
        children = node.get("children", [])

        # If no children or at a mol node with no heterocycles, we're at the end of a potential path
        if len(children) == 0 and node["type"] == "mol" and len(current_path) > 0:
            # Check if this is a starting material
            if node.get("in_stock", False):
                # Save this path
                paths.append(list(current_path))
                print(f"Completed path with {len(current_path)} heterocycle-containing molecules")

        # Continue traversal
        path_length_before = len(current_path)
        for child in children:
            dfs_traverse(child, depth + 1)

        # Backtrack: remove any molecules added at this level
        while len(current_path) > 0 and current_path[-1]["depth"] >= depth:
            current_path.pop()

    # Start traversal
    dfs_traverse(route)

    # Analyze paths to find linear heterocycle synthesis strategy
    print(f"Found {len(paths)} potential heterocycle synthesis paths")

    for i, path in enumerate(paths):
        print(f"Path {i+1}: {len(path)} heterocycle-containing molecules")

        # Check if this path has sequential heterocycle transformations
        if len(path) >= 3:
            # Sort by depth to ensure we're following the synthetic direction
            sorted_path = sorted(path, key=lambda x: x["depth"])

            # Check for sequential transformations (at least 3)
            sequential_transformations = 0
            for j in range(len(sorted_path) - 1):
                current_mol = sorted_path[j]
                next_mol = sorted_path[j + 1]

                # If depths differ by 2, there's a reaction node between them
                if next_mol["depth"] - current_mol["depth"] == 2:
                    sequential_transformations += 1
                    print(
                        f"  Sequential transformation {sequential_transformations}: "
                        f"{current_mol['heterocycles']} → {next_mol['heterocycles']}"
                    )

            if (
                sequential_transformations >= 2
            ):  # At least 3 molecules with 2 transformations between them
                print(
                    f"Found linear heterocycle synthesis strategy with {sequential_transformations} sequential transformations"
                )
                return True

    print("No linear heterocycle synthesis strategy found")
    return False
