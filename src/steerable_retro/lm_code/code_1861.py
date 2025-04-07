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
    This function detects a strategy where a heterocycle (specifically benzimidazole)
    is formed in the middle of the synthesis, followed by late-stage cross-coupling.
    """
    # Track if we found the key features and their depths
    heterocycle_formation_depth = None
    cross_coupling_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_depth, cross_coupling_depth

        if node["type"] == "reaction":
            # Extract reaction information
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for heterocycle formation
            # Check for benzimidazole formation
            if checker.check_ring("benzimidazole", product_smiles):
                reactants_have_benzimidazole = any(
                    checker.check_ring("benzimidazole", r) for r in reactants_smiles
                )
                if not reactants_have_benzimidazole:
                    heterocycle_formation_depth = depth
                    print(f"Found benzimidazole formation at depth {depth}")

            # Check for other relevant heterocycles
            for ring in [
                "benzoxazole",
                "benzothiazole",
                "imidazole",
                "oxazole",
                "thiazole",
                "triazole",
            ]:
                if checker.check_ring(ring, product_smiles):
                    reactants_have_ring = any(checker.check_ring(ring, r) for r in reactants_smiles)
                    if not reactants_have_ring:
                        heterocycle_formation_depth = depth
                        print(f"Found {ring} formation at depth {depth}")

            # Check for cross-coupling reactions
            # Check for Suzuki coupling
            if (
                checker.check_reaction("Suzuki coupling with boronic acids", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic acids OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with boronic esters OTf", rsmi)
                or checker.check_reaction("Suzuki coupling with sulfonic esters", rsmi)
            ):
                cross_coupling_depth = depth
                print(f"Found Suzuki cross-coupling at depth {depth}")

            # Check for other cross-coupling reactions
            elif (
                checker.check_reaction("Negishi coupling", rsmi)
                or checker.check_reaction("Stille reaction_aryl", rsmi)
                or checker.check_reaction("Stille reaction_vinyl", rsmi)
                or checker.check_reaction("Hiyama-Denmark Coupling", rsmi)
                or checker.check_reaction("Kumada cross-coupling", rsmi)
                or checker.check_reaction("Sonogashira alkyne_aryl halide", rsmi)
                or checker.check_reaction("Sonogashira acetylene_aryl halide", rsmi)
            ):
                cross_coupling_depth = depth
                print(f"Found other cross-coupling at depth {depth}")

            # Fallback check using functional groups if reaction check fails
            elif depth <= 3:  # Only check for late-stage reactions
                has_boronic = any(
                    checker.check_fg("Boronic acid", r) or checker.check_fg("Boronic ester", r)
                    for r in reactants_smiles
                )
                has_aryl_halide = any(
                    checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                )
                has_metal = any(
                    checker.check_fg("Zinc halide", r) or checker.check_fg("Magnesium halide", r)
                    for r in reactants_smiles
                )

                if (has_boronic or has_metal) and has_aryl_halide:
                    cross_coupling_depth = depth
                    print(f"Found cross-coupling (detected by functional groups) at depth {depth}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if both key features were found and in the correct order
    # (heterocycle formation should happen at a greater depth than cross-coupling)
    if heterocycle_formation_depth is not None and cross_coupling_depth is not None:
        print(
            f"Heterocycle formation depth: {heterocycle_formation_depth}, Cross-coupling depth: {cross_coupling_depth}"
        )
        return heterocycle_formation_depth > cross_coupling_depth

    return False
