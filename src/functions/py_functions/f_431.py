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
    This function detects a synthetic strategy where a sulfonamide formation occurs
    in the late stage of the synthesis (low depth).
    """
    sulfonamide_formed = False
    sulfonamide_depth = float("inf")

    def dfs_traverse(node, depth=0):
        nonlocal sulfonamide_formed, sulfonamide_depth

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]

                # Check if this is a sulfonamide formation reaction using the checker function
                is_primary_sulfonamide = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                )
                is_secondary_sulfonamide = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )

                if is_primary_sulfonamide or is_secondary_sulfonamide:
                    print(f"Found sulfonamide formation at depth {depth}, rsmi: {rsmi}")
                    sulfonamide_formed = True
                    sulfonamide_depth = min(sulfonamide_depth, depth)

                # Fallback check if the specific reaction types aren't detected
                if not (is_primary_sulfonamide or is_secondary_sulfonamide):
                    reactants = rsmi.split(">")[0].split(".")
                    product = rsmi.split(">")[-1]

                    # Check for sulfonyl chloride in reactants
                    reactant_has_sulfonyl_chloride = any(
                        checker.check_fg("Sulfonyl halide", r) for r in reactants
                    )

                    # Check for sulfonamide in product
                    product_has_sulfonamide = checker.check_fg("Sulfonamide", product)

                    if reactant_has_sulfonyl_chloride and product_has_sulfonamide:
                        print(
                            f"Found sulfonamide formation via FG check at depth {depth}, rsmi: {rsmi}"
                        )
                        sulfonamide_formed = True
                        sulfonamide_depth = min(sulfonamide_depth, depth)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Consider it late-stage if it's in the first few steps of the synthesis (depth <= 2)
    is_late_stage = sulfonamide_formed and sulfonamide_depth <= 2
    print(
        f"Sulfonamide formed: {sulfonamide_formed}, at depth: {sulfonamide_depth}, is late stage: {is_late_stage}"
    )
    return is_late_stage
