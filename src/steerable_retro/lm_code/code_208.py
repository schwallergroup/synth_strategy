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
    This function detects a synthetic strategy involving late-stage reductive amination
    to form a C-N bond as the final step in the synthesis.
    """
    # Track if we found reductive amination at depth 0 or 1 (final steps)
    found_late_stage_reductive_amination = False

    def dfs_traverse(node, current_depth=0):
        nonlocal found_late_stage_reductive_amination

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                # Handle depth properly - use current_depth if metadata depth is not available
                depth = current_depth
                if "depth" in node.get("metadata", {}):
                    depth = int(node["metadata"]["depth"])

                # Check if this is a late-stage step (depth 0, 1, or -1)
                if depth <= 1 or depth == -1:
                    print(f"Examining reaction at depth {depth}: {rsmi}")

                    # Check if this is a reductive amination reaction
                    is_reductive_amination_aldehyde = checker.check_reaction(
                        "Reductive amination with aldehyde", rsmi
                    )
                    is_reductive_amination_ketone = checker.check_reaction(
                        "Reductive amination with ketone", rsmi
                    )
                    is_reductive_amination_alcohol = checker.check_reaction(
                        "Reductive amination with alcohol", rsmi
                    )

                    if (
                        is_reductive_amination_aldehyde
                        or is_reductive_amination_ketone
                        or is_reductive_amination_alcohol
                    ):
                        print(
                            f"Detected reductive amination reaction: aldehyde={is_reductive_amination_aldehyde}, ketone={is_reductive_amination_ketone}, alcohol={is_reductive_amination_alcohol}"
                        )

                        # Get reactants and product
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check for aldehyde/ketone/alcohol and amine in reactants
                        has_aldehyde = False
                        has_ketone = False
                        has_alcohol = False
                        has_primary_amine = False
                        has_secondary_amine = False

                        for reactant in reactants:
                            if checker.check_fg("Aldehyde", reactant):
                                has_aldehyde = True
                                print(f"Found aldehyde in reactant: {reactant}")

                            if checker.check_fg("Ketone", reactant):
                                has_ketone = True
                                print(f"Found ketone in reactant: {reactant}")

                            if checker.check_fg("Primary alcohol", reactant) or checker.check_fg(
                                "Secondary alcohol", reactant
                            ):
                                has_alcohol = True
                                print(f"Found alcohol in reactant: {reactant}")

                            if checker.check_fg("Primary amine", reactant):
                                has_primary_amine = True
                                print(f"Found primary amine in reactant: {reactant}")

                            if checker.check_fg("Secondary amine", reactant):
                                has_secondary_amine = True
                                print(f"Found secondary amine in reactant: {reactant}")

                        # Check for secondary or tertiary amine in product
                        has_secondary_amine_product = checker.check_fg("Secondary amine", product)
                        has_tertiary_amine_product = checker.check_fg("Tertiary amine", product)

                        if has_secondary_amine_product:
                            print(f"Found secondary amine in product: {product}")

                        if has_tertiary_amine_product:
                            print(f"Found tertiary amine in product: {product}")

                        # Confirm we have the right reactants and products for reductive amination
                        primary_to_secondary = has_primary_amine and has_secondary_amine_product
                        secondary_to_tertiary = has_secondary_amine and has_tertiary_amine_product

                        if (has_aldehyde or has_ketone or has_alcohol) and (
                            primary_to_secondary or secondary_to_tertiary
                        ):
                            found_late_stage_reductive_amination = True
                            print(
                                f"Confirmed late-stage reductive amination strategy at depth {depth}"
                            )

                    # Additional check for reductive amination pattern if the reaction checker didn't catch it
                    elif not found_late_stage_reductive_amination:
                        # Get reactants and product
                        reactants = rsmi.split(">")[0].split(".")
                        product = rsmi.split(">")[-1]

                        # Check for aldehyde/ketone and amine in reactants
                        has_aldehyde = any(checker.check_fg("Aldehyde", r) for r in reactants)
                        has_ketone = any(checker.check_fg("Ketone", r) for r in reactants)
                        has_primary_amine = any(
                            checker.check_fg("Primary amine", r) for r in reactants
                        )
                        has_secondary_amine = any(
                            checker.check_fg("Secondary amine", r) for r in reactants
                        )

                        # Check for secondary or tertiary amine in product
                        has_secondary_amine_product = checker.check_fg("Secondary amine", product)
                        has_tertiary_amine_product = checker.check_fg("Tertiary amine", product)

                        # Look for pattern where aldehyde/ketone + amine → amine product
                        if (has_aldehyde or has_ketone) and (
                            (has_primary_amine and has_secondary_amine_product)
                            or (has_secondary_amine and has_tertiary_amine_product)
                        ):

                            # Check if the product has a new C-N bond that wasn't in the reactants
                            # This is a simplified check - in a real implementation, you'd need to analyze
                            # the atom mapping more carefully
                            found_late_stage_reductive_amination = True
                            print(f"Detected reductive amination pattern at depth {depth}")

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, current_depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Late-stage reductive amination strategy detected: {found_late_stage_reductive_amination}"
    )

    return found_late_stage_reductive_amination
