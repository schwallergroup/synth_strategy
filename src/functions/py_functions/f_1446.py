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
    This function detects a strategy involving late-stage nitration
    (nitro group introduction in the final steps of synthesis).
    """
    nitration_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal nitration_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if this is a nitration reaction
                is_nitration = (
                    checker.check_reaction("Aromatic nitration with HNO3", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO3 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with NO2 salt", rsmi)
                    or checker.check_reaction("Aromatic nitration with alkyl NO2", rsmi)
                    or checker.check_reaction("Non-aromatic nitration with HNO3", rsmi)
                )

                # If it's a nitration reaction, verify nitro group is introduced
                if is_nitration:
                    # Verify nitro group is in product but not in reactants
                    product_has_nitro = checker.check_fg("Nitro group", product)
                    reactants_have_nitro = any(
                        checker.check_fg("Nitro group", r) for r in reactants
                    )

                    if product_has_nitro and not reactants_have_nitro:
                        print(f"Detected nitration at depth {depth}, rsmi: {rsmi}")
                        nitration_depth = depth

                # Alternative detection method if reaction check fails
                if nitration_depth is None:
                    product_has_nitro = checker.check_fg("Nitro group", product)
                    reactants_have_nitro = any(
                        checker.check_fg("Nitro group", r) for r in reactants
                    )

                    if product_has_nitro and not reactants_have_nitro:
                        print(
                            f"Detected potential nitration at depth {depth} by FG analysis, rsmi: {rsmi}"
                        )
                        nitration_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if strategy is present (nitration in first 25% of steps)
    # In retrosynthetic analysis, lower depth means later stage in synthesis
    strategy_present = (
        nitration_depth is not None and nitration_depth <= max_depth * 0.25
    )

    print(f"Nitration depth: {nitration_depth}")
    print(f"Max depth: {max_depth}")
    print(f"Strategy detected: {strategy_present}")

    return strategy_present
