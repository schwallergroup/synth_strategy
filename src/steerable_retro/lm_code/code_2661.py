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
    This function detects a pattern of C-S bond formation followed by N-S bond formation.
    Specifically, it looks for sulfonyl chloride formation followed by sulfonamide formation.
    """
    # Track reactions and their depths
    sulfonyl_chloride_formation = None  # Will store (depth, reaction_smiles)
    sulfonamide_formation = None  # Will store (depth, reaction_smiles)

    def dfs_traverse(node, depth=0):
        nonlocal sulfonyl_chloride_formation, sulfonamide_formation

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for sulfonyl chloride formation (C-S bond)
                if checker.check_reaction("Aromatic sulfonyl chlorination", rsmi):
                    # Verify product has sulfonyl chloride but reactants don't
                    if checker.check_fg("Sulfonyl halide", product):
                        has_sulfonyl_in_reactants = any(
                            checker.check_fg("Sulfonyl halide", r) for r in reactants
                        )
                        if not has_sulfonyl_in_reactants:
                            print(f"Detected sulfonyl chloride formation at depth {depth}: {rsmi}")
                            sulfonyl_chloride_formation = (depth, rsmi)

                # Check for sulfonamide formation (N-S bond)
                if checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                ):
                    # Verify product has sulfonamide but reactants don't
                    if checker.check_fg("Sulfonamide", product):
                        has_sulfonamide_in_reactants = any(
                            checker.check_fg("Sulfonamide", r) for r in reactants
                        )
                        if not has_sulfonamide_in_reactants:
                            print(f"Detected sulfonamide formation at depth {depth}: {rsmi}")
                            sulfonamide_formation = (depth, rsmi)

                # Alternative check for C-S bond formation if specific reaction check fails
                if sulfonyl_chloride_formation is None:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and checker.check_fg("Sulfonyl halide", product):
                        # Check if reactants don't have sulfonyl halide
                        has_sulfonyl_in_reactants = any(
                            checker.check_fg("Sulfonyl halide", r) for r in reactants
                        )
                        if not has_sulfonyl_in_reactants:
                            print(
                                f"Detected sulfonyl halide formation (alternative) at depth {depth}: {rsmi}"
                            )
                            sulfonyl_chloride_formation = (depth, rsmi)

                # Alternative check for N-S bond formation if specific reaction check fails
                if sulfonamide_formation is None:
                    product_mol = Chem.MolFromSmiles(product)
                    if product_mol and checker.check_fg("Sulfonamide", product):
                        # Check if reactants don't have sulfonamide
                        has_sulfonamide_in_reactants = any(
                            checker.check_fg("Sulfonamide", r) for r in reactants
                        )
                        if not has_sulfonamide_in_reactants:
                            print(
                                f"Detected sulfonamide formation (alternative) at depth {depth}: {rsmi}"
                            )
                            sulfonamide_formation = (depth, rsmi)
            except Exception as e:
                print(f"Error in processing molecules for sulfonylation pattern check: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Call dfs_traverse on the root node
    dfs_traverse(route)

    # Check if both reactions were found and in the correct order
    result = False
    if sulfonyl_chloride_formation and sulfonamide_formation:
        # Ensure C-S formation happens before N-S formation (lower depth means later in synthesis)
        if sulfonyl_chloride_formation[0] > sulfonamide_formation[0]:
            result = True

    print(f"Sulfonylation pattern (C-S followed by N-S): {result}")
    print(
        f"  - Sulfonyl chloride formation: {sulfonyl_chloride_formation[1] if sulfonyl_chloride_formation else None}"
    )
    print(
        f"  - Sulfonamide formation: {sulfonamide_formation[1] if sulfonamide_formation else None}"
    )

    return result
