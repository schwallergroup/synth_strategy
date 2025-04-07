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
    This function detects if a synthetic route maintains a nitrile functional group
    through multiple steps before transforming it in a late stage.
    """
    # Track nitrile presence and transformations
    nitrile_presence = []
    nitrile_transformations = []
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_nitrile = checker.check_fg("Nitrile", mol_smiles)

            # Store information about nitrile presence
            nitrile_presence.append((depth, has_nitrile, mol_smiles))

            # Process children (reactions)
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

        elif node["type"] == "reaction":
            # Check if this reaction transforms a nitrile
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if any reactant has nitrile but product doesn't
                reactants_with_nitrile = any(
                    checker.check_fg("Nitrile", r) for r in reactants_smiles
                )
                product_has_nitrile = checker.check_fg("Nitrile", product_smiles)

                if reactants_with_nitrile and not product_has_nitrile:
                    # This reaction transforms a nitrile
                    nitrile_transformations.append((depth, rsmi))
                    print(f"Nitrile transformation at depth {depth}: {rsmi}")
            except Exception as e:
                print(f"Error analyzing reaction: {e}")

            # Process children (reactants)
            for child in node.get("children", []):
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Sort by depth to analyze the synthesis route
    nitrile_presence.sort(key=lambda x: x[0])
    nitrile_transformations.sort(key=lambda x: x[0])

    # Check if nitrile is present in early stages (higher depth)
    early_stages = [item for item in nitrile_presence if item[0] > max_depth // 2]
    nitrile_in_early_stages = any(has_nitrile for _, has_nitrile, _ in early_stages)

    # Check if nitrile is maintained through multiple steps
    depths_with_nitrile = set(
        depth for depth, has_nitrile, _ in nitrile_presence if has_nitrile
    )
    consecutive_depths = len(depths_with_nitrile) > 1

    # Check if nitrile is transformed in late stage (lower depth)
    late_stage_transformation = any(
        depth <= max_depth // 2 for depth, _ in nitrile_transformations
    )

    print(f"Nitrile in early stages: {nitrile_in_early_stages}")
    print(f"Nitrile maintained through steps: {consecutive_depths}")
    print(f"Late stage transformation: {late_stage_transformation}")
    print(f"Max depth: {max_depth}")

    # Route maintains nitrile until late stage if:
    # 1. Nitrile is present in early stages
    # 2. Nitrile is maintained through multiple steps
    # 3. Nitrile is transformed in a late stage
    result = (
        nitrile_in_early_stages and consecutive_depths and late_stage_transformation
    )

    if result:
        print("Route maintains nitrile until late stage transformation")
    else:
        print("Route does not maintain nitrile until late stage")

    return result
