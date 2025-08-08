#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

root_data = "/home/dparm/steerable_retro/data"

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
    This function detects a synthetic strategy involving late-stage sulfonamide formation
    with an isoxazole fragment.
    """
    found_pattern = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern

        if found_pattern:
            return

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this is a late-stage reaction (depth 0 or 1)
            if depth <= 1:
                print(f"Checking late-stage reaction at depth {depth}: {rsmi}")

                # Check if this is a sulfonamide formation reaction
                is_sulfonamide_reaction = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                ) or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )

                # If not detected by reaction checker, try manual detection
                if not is_sulfonamide_reaction:
                    print("Reaction not detected by reaction checker, trying manual detection...")

                    # Check for sulfonyl halide in reactants
                    sulfonyl_halide_present = False
                    for reactant in reactants:
                        if checker.check_fg("Sulfonyl halide", reactant):
                            sulfonyl_halide_present = True
                            print(f"Found sulfonyl halide in reactant: {reactant}")
                            break

                    # Check for amine in reactants
                    amine_present = False
                    for reactant in reactants:
                        if checker.check_fg("Primary amine", reactant) or checker.check_fg(
                            "Secondary amine", reactant
                        ):
                            amine_present = True
                            print(f"Found amine in reactant: {reactant}")
                            break

                    # Check for sulfonamide in product but not in reactants
                    sulfonamide_in_product = checker.check_fg("Sulfonamide", product)
                    sulfonamide_in_reactants = any(
                        checker.check_fg("Sulfonamide", r) for r in reactants
                    )

                    if (
                        sulfonyl_halide_present
                        and amine_present
                        and sulfonamide_in_product
                        and not sulfonamide_in_reactants
                    ):
                        print("Manually detected sulfonamide formation reaction")
                        is_sulfonamide_reaction = True

                if is_sulfonamide_reaction:
                    print("Found sulfonamide formation reaction")

                    # Check if isoxazole is present in the product
                    if checker.check_ring("isoxazole", product):
                        print("Found isoxazole in product")

                        # Verify isoxazole is present in at least one reactant
                        isoxazole_in_reactants = False
                        for reactant in reactants:
                            if checker.check_ring("isoxazole", reactant):
                                isoxazole_in_reactants = True
                                print(f"Found isoxazole in reactant: {reactant}")
                                break

                        if not isoxazole_in_reactants:
                            print("No isoxazole found in reactants, skipping")
                            return

                        # Check for sulfonyl chloride in reactants
                        sulfonyl_chloride_present = False
                        for reactant in reactants:
                            if checker.check_fg("Sulfonyl halide", reactant):
                                sulfonyl_chloride_present = True
                                print(f"Found sulfonyl chloride in reactant: {reactant}")
                                break

                        if not sulfonyl_chloride_present:
                            print("No sulfonyl chloride found in reactants, skipping")
                            return

                        # Verify sulfonamide was actually formed (not present in reactants)
                        sulfonamide_in_reactants = False
                        for reactant in reactants:
                            if checker.check_fg("Sulfonamide", reactant):
                                sulfonamide_in_reactants = True
                                print(f"Found sulfonamide already in reactant: {reactant}")
                                break

                        if sulfonamide_in_reactants:
                            print("Sulfonamide already present in reactants, skipping")
                            return

                        if checker.check_fg("Sulfonamide", product):
                            print("Confirmed sulfonamide formation with isoxazole")
                            found_pattern = True
                        else:
                            print("No sulfonamide found in product, skipping")
                    else:
                        print("No isoxazole found in product, skipping")
                else:
                    print("Not a sulfonamide formation reaction, skipping")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return found_pattern
