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
    This function detects a strategy involving late-stage coupling of heterocycles via hydrazone formation.
    It checks for:
    1. Hydrazone formation in the final steps
    2. Connection of two heterocyclic systems (quinoline and pyridine)
    """
    final_coupling_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal final_coupling_detected

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            # Check reactions in the late stage (depth 0-2)
            if depth <= 2:  # Late stage reactions
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Examining reaction at depth {depth}: {rsmi}")

                # Check if this is a hydrazone formation reaction
                is_hydrazone_formation = checker.check_reaction(
                    "Ketone/aldehyde to hydrazone", rsmi
                )

                # Verify product contains hydrazone
                has_hydrazone_product = checker.check_fg("Hydrazone", product_smiles)

                print(f"Is hydrazone formation: {is_hydrazone_formation}")
                print(f"Product has hydrazone: {has_hydrazone_product}")

                # Check if any reactant already has hydrazone
                reactants_with_hydrazone = [
                    r for r in reactants_smiles if checker.check_fg("Hydrazone", r)
                ]

                if (
                    is_hydrazone_formation or has_hydrazone_product
                ) and not reactants_with_hydrazone:
                    print(f"Found potential hydrazone formation at depth {depth}")

                    # Check if reactants contain heterocyclic systems
                    quinoline_reactants = []
                    pyridine_reactants = []
                    hydrazine_reactants = []
                    carbonyl_reactants = []

                    for reactant in reactants_smiles:
                        # Check for heterocycles
                        if checker.check_ring("quinoline", reactant):
                            quinoline_reactants.append(reactant)
                            print(f"Found quinoline in reactant: {reactant}")

                        if checker.check_ring("pyridine", reactant):
                            pyridine_reactants.append(reactant)
                            print(f"Found pyridine in reactant: {reactant}")

                        # Check for functional groups involved in hydrazone formation
                        if checker.check_fg("Hydrazine", reactant) or checker.check_fg(
                            "Acylhydrazine", reactant
                        ):
                            hydrazine_reactants.append(reactant)
                            print(f"Found hydrazine/hydrazide in reactant: {reactant}")

                        if checker.check_fg("Aldehyde", reactant) or checker.check_fg(
                            "Ketone", reactant
                        ):
                            carbonyl_reactants.append(reactant)
                            print(f"Found carbonyl in reactant: {reactant}")

                    # Filter out quinoline reactants from pyridine list to avoid double counting
                    pyridine_only_reactants = [
                        r for r in pyridine_reactants if r not in quinoline_reactants
                    ]

                    has_quinoline = len(quinoline_reactants) > 0
                    has_pyridine = (
                        len(pyridine_only_reactants) > 0 or len(pyridine_reactants) > 0
                    )
                    has_hydrazine = len(hydrazine_reactants) > 0
                    has_carbonyl = len(carbonyl_reactants) > 0

                    print(
                        f"Has quinoline: {has_quinoline}, Has pyridine: {has_pyridine}"
                    )
                    print(
                        f"Has hydrazine: {has_hydrazine}, Has carbonyl: {has_carbonyl}"
                    )

                    # Verify the heterocycles are connected via the hydrazone in the product
                    if (
                        has_hydrazone_product
                        and (has_quinoline or has_pyridine)
                        and (has_hydrazine or has_carbonyl)
                    ):
                        print("Found heterocycle with required functional groups")

                        # Check if the product contains both heterocycles
                        product_has_quinoline = checker.check_ring(
                            "quinoline", product_smiles
                        )
                        product_has_pyridine = checker.check_ring(
                            "pyridine", product_smiles
                        )

                        print(f"Product has quinoline: {product_has_quinoline}")
                        print(f"Product has pyridine: {product_has_pyridine}")

                        if product_has_quinoline and product_has_pyridine:
                            print("Product contains both heterocycles")

                            # Verify the product has a hydrazone that wasn't in any reactant
                            if not any(
                                checker.check_fg("Hydrazone", r)
                                for r in reactants_smiles
                            ):
                                final_coupling_detected = True
                                print(
                                    f"Detected late-stage hydrazone coupling of heterocycles at depth {depth}"
                                )

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return final_coupling_detected
