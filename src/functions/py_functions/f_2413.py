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
    This function detects oxazole ring formation in the early stages of synthesis.
    """
    early_stage_oxazole_found = False
    max_depth = 0

    # First pass to determine the maximum depth
    def get_max_depth(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)
        for child in node.get("children", []):
            get_max_depth(child, depth + 1)

    get_max_depth(route)
    print(f"Maximum depth of synthesis route: {max_depth}")

    # Second pass to find oxazole formation in early stage
    def dfs_traverse(node, depth=0):
        nonlocal early_stage_oxazole_found

        if node["type"] == "reaction":
            # Extract reactants and product
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product contains oxazole
                product_has_oxazole = checker.check_ring("oxazole", product_smiles)

                # Check if any reactant contains oxazole
                reactants_have_oxazole = any(
                    checker.check_ring("oxazole", reactant)
                    for reactant in reactants_smiles
                )

                # Check if this is an oxazole formation reaction
                if product_has_oxazole and not reactants_have_oxazole:
                    print(f"Oxazole formation detected at depth {depth}")

                    # Check if this is early stage (depth > max_depth/2)
                    if depth > (max_depth / 2):
                        print(
                            f"Early stage oxazole formation confirmed at depth {depth} (max depth: {max_depth})"
                        )
                        early_stage_oxazole_found = True

                # Also check for specific oxazole formation reactions
                if (
                    checker.check_reaction("benzoxazole formation from aldehyde", rsmi)
                    or checker.check_reaction(
                        "benzoxazole formation from acyl halide", rsmi
                    )
                    or checker.check_reaction(
                        "benzoxazole formation from ester/carboxylic acid", rsmi
                    )
                    or checker.check_reaction(
                        "benzoxazole formation (intramolecular)", rsmi
                    )
                ):
                    print(
                        f"Specific oxazole formation reaction detected at depth {depth}"
                    )
                    if depth > (max_depth / 2):
                        print(
                            f"Early stage oxazole formation reaction confirmed at depth {depth}"
                        )
                        early_stage_oxazole_found = True

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    print(f"Early stage oxazole formation found: {early_stage_oxazole_found}")
    return early_stage_oxazole_found
