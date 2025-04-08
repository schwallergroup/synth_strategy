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
    This function detects if the final product contains multiple nitrogen heterocycles.
    """
    has_multiple_heterocycles = False

    def dfs_traverse(node, depth=0):
        nonlocal has_multiple_heterocycles

        # Check only the root node (final product)
        if depth == 0 and node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            print(f"Analyzing final product: {smiles}")

            # List of nitrogen-containing heterocycles to check
            nitrogen_heterocycles = [
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
                "carbazole",
                "acridine",
                "benzimidazole",
                "indazole",
                "benzotriazole",
            ]

            # Count different heterocycle types
            heterocycle_types = 0
            found_rings = []

            for ring in nitrogen_heterocycles:
                if checker.check_ring(ring, smiles):
                    heterocycle_types += 1
                    found_rings.append(ring)
                    print(f"{ring.capitalize()} ring detected in final product")

            print(f"Total nitrogen heterocycle types found: {heterocycle_types}")
            print(f"Found rings: {', '.join(found_rings)}")

            if heterocycle_types >= 2:
                has_multiple_heterocycles = True
                print("Multiple nitrogen heterocycle types detected in final product")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_multiple_heterocycles
