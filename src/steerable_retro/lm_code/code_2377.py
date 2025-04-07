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
    Detects if the route contains protection of an aldehyde as a cyclic acetal.
    """
    has_acetal_protection = False

    def dfs_traverse(node):
        nonlocal has_acetal_protection

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            print(f"Examining reaction: {rsmi}")

            # Check if this is an acetalization reaction (either perspective)
            is_acetalization = checker.check_reaction(
                "Aldehyde or ketone acetalization", rsmi
            ) or checker.check_reaction("Diol acetalization", rsmi)

            if is_acetalization:
                print(f"Found acetalization reaction: {rsmi}")

                # Verify that an aldehyde is being protected (not a ketone)
                aldehyde_reactant = None
                for reactant in reactants:
                    if checker.check_fg("Aldehyde", reactant):
                        print(f"Found aldehyde in reactants: {reactant}")
                        aldehyde_reactant = reactant
                        break

                if aldehyde_reactant:
                    # Check if the product contains an acetal/ketal
                    if checker.check_fg("Acetal/Ketal", product):
                        print(f"Found acetal in product: {product}")

                        # Verify it's a cyclic acetal by checking for dioxolane or dioxane rings
                        if (
                            checker.check_ring("dioxolane", product)
                            or checker.check_ring("dioxane", product)
                            or checker.check_ring("dioxepane", product)
                        ):
                            print("Confirmed cyclic acetal structure")

                            # Verify the aldehyde is consumed in the reaction
                            # In the forward direction, the aldehyde should disappear
                            # and be replaced by an acetal
                            if not checker.check_fg("Aldehyde", product):
                                print("Confirmed aldehyde is consumed in the reaction")
                                has_acetal_protection = True

            # Also check for the reverse reaction (deprotection) which in retrosynthesis would be protection
            is_acetal_hydrolysis = checker.check_reaction("Acetal hydrolysis to aldehyde", rsmi)

            if is_acetal_hydrolysis:
                print(f"Found acetal hydrolysis reaction (reverse of protection): {rsmi}")

                # In deprotection, the product should have an aldehyde
                if checker.check_fg("Aldehyde", product):
                    print(f"Found aldehyde in product (deprotection): {product}")

                    # And at least one reactant should have a cyclic acetal
                    for reactant in reactants:
                        print(f"Checking reactant for acetal: {reactant}")
                        if checker.check_fg("Acetal/Ketal", reactant):
                            print(f"Found acetal in reactant: {reactant}")

                            # Verify it's a cyclic acetal
                            if (
                                checker.check_ring("dioxolane", reactant)
                                or checker.check_ring("dioxane", reactant)
                                or checker.check_ring("dioxepane", reactant)
                            ):
                                print("Confirmed cyclic acetal structure in deprotection")
                                has_acetal_protection = True
                                break

            # Special case: Check for a reaction that looks like acetal hydrolysis
            # by examining the structural changes directly
            if not has_acetal_protection:
                # Check if product contains an aldehyde
                if checker.check_fg("Aldehyde", product):
                    print(f"Found aldehyde in product: {product}")

                    # Check if any reactant contains a cyclic acetal structure
                    for reactant in reactants:
                        if (
                            "C1C[O" in reactant or "C1CCO" in reactant
                        ):  # Simple pattern for dioxolane
                            print(f"Found potential cyclic acetal in reactant: {reactant}")

                            # Verify with more specific checks
                            if (
                                checker.check_ring("dioxolane", reactant)
                                or checker.check_ring("dioxane", reactant)
                                or checker.check_ring("dioxepane", reactant)
                            ):
                                print("Confirmed cyclic acetal structure")

                                # Check if the acetal carbon is connected to the aromatic system
                                # This is a heuristic check for the acetal protecting an aldehyde
                                if "[CH" in reactant and "[c" in reactant:
                                    print("Detected potential acetal protecting an aldehyde")
                                    has_acetal_protection = True
                                    break

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Final result: {has_acetal_protection}")
    return has_acetal_protection
