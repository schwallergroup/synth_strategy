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
    This function detects a strategy with early cyclopropane formation and
    late-stage nitrile hydrolysis to carboxylic acid.
    """
    cyclopropane_formed_early = False
    nitrile_hydrolysis_late = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal cyclopropane_formed_early, nitrile_hydrolysis_late, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for cyclopropane formation
            if depth >= 2:  # Early in synthesis (high depth)
                # Check if cyclopropane was formed in this step
                if checker.check_ring("cyclopropane", product):
                    if not any(checker.check_ring("cyclopropane", r) for r in reactants):
                        print(f"Cyclopropane formation detected at depth {depth}")
                        cyclopropane_formed_early = True

            # Check for nitrile hydrolysis to carboxylic acid
            if depth <= 1:  # Late stage (final or penultimate step)
                # Check for nitrile hydrolysis reaction
                if checker.check_reaction("Oxidation of nitrile to carboxylic acid", rsmi):
                    print(f"Nitrile hydrolysis to carboxylic acid detected at depth {depth}")
                    nitrile_hydrolysis_late = True
                else:
                    # Fallback check if reaction type not directly matched
                    if checker.check_fg("Carboxylic acid", product):
                        if any(checker.check_fg("Nitrile", r) for r in reactants):
                            print(
                                f"Nitrile hydrolysis to carboxylic acid detected at depth {depth} (via FG check)"
                            )
                            nitrile_hydrolysis_late = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Cyclopropane formed early: {cyclopropane_formed_early}")
    print(f"Nitrile hydrolysis late: {nitrile_hydrolysis_late}")
    print(f"Max depth: {max_depth}")

    return cyclopropane_formed_early and nitrile_hydrolysis_late
