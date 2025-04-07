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
    This function detects if a heterocycle is formed
    in the middle stages of the synthesis (not at the beginning or end).
    """
    # Track if we found heterocycle formation
    heterocycle_formed = False
    total_depth = 0
    max_depth = 0
    heterocycle_depth = None

    # List of heterocycles to check for
    heterocycles = [
        "isoxazole",
        "oxazole",
        "thiazole",
        "pyrrole",
        "furan",
        "pyridine",
        "imidazole",
        "triazole",
        "tetrazole",
        "pyrazole",
        "oxadiazole",
        "thiadiazole",
        "pyrimidine",
        "pyrazine",
        "pyridazine",
    ]

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formed, total_depth, max_depth, heterocycle_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check which heterocycles are formed in this reaction
            formed_heterocycles = []
            for ring in heterocycles:
                # Check if ring is in product but not in any reactant
                if checker.check_ring(ring, product_smiles):
                    if not any(checker.check_ring(ring, r) for r in reactants_smiles):
                        formed_heterocycles.append(ring)

            # If any heterocycle is formed in this step
            if formed_heterocycles:
                heterocycle_formed = True
                heterocycle_depth = depth
                print(
                    f"Heterocycle formation detected at depth {depth}: {', '.join(formed_heterocycles)}"
                )

                # Also check if this is a known heterocycle formation reaction
                if checker.check_reaction("Formation of NOS Heterocycles", rsmi):
                    print(f"Confirmed as NOS Heterocycle formation reaction")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    total_depth = max_depth

    # Check if heterocycle formation occurred in the middle stages
    # (not at the beginning or end of the synthesis)
    if heterocycle_formed and heterocycle_depth is not None:
        # Consider middle stage as between 25% and 75% of the total synthesis depth
        lower_bound = total_depth * 0.25
        upper_bound = total_depth * 0.75

        is_mid_stage = lower_bound <= heterocycle_depth <= upper_bound
        print(f"Heterocycle formed at depth {heterocycle_depth} out of total depth {total_depth}")
        print(f"Middle stage bounds: {lower_bound} to {upper_bound}")
        print(f"Is mid-stage heterocycle formation: {is_mid_stage}")

        return is_mid_stage

    return False
