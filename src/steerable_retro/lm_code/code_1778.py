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
    This function detects a 'redox dance' pattern where a sequence of
    oxidation-reduction-oxidation steps are used in the synthesis.

    In retrosynthesis, we traverse from the target molecule backward, so the depths
    increase as we go earlier in the synthesis. The redox dance pattern should be
    detected as: oxidation (late stage) → reduction (middle stage) → oxidation (early stage)
    """
    # Store reaction steps with their depths
    oxidation_steps = []
    reduction_steps = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # In retrosynthesis, the product is what we start with and reactants are what we're making
                # So for oxidation: product (alcohol) -> reactant (carbonyl)
                # For reduction: product (carbonyl) -> reactant (alcohol)

                # Check for oxidation reactions (alcohol to carbonyl, etc.)
                if (
                    checker.check_reaction("Oxidation of aldehydes to carboxylic acids", rsmi)
                    or checker.check_reaction("Oxidation of ketone to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of alcohol to carboxylic acid", rsmi)
                    or checker.check_reaction(
                        "Oxidation or Dehydrogenation of Alcohols to Aldehydes and Ketones", rsmi
                    )
                    or checker.check_reaction("Oxidation of alkene to carboxylic acid", rsmi)
                    or checker.check_reaction("Oxidation of alkene to aldehyde", rsmi)
                    or checker.check_reaction("Oxidation of boronic acids", rsmi)
                    or checker.check_reaction("Oxidation of boronic esters", rsmi)
                ):
                    oxidation_steps.append(depth)
                    print(f"Detected oxidation at depth {depth}: {rsmi}")

                # Check for reduction reactions (carbonyl to alcohol, etc.)
                if (
                    checker.check_reaction("Reduction of aldehydes and ketones to alcohols", rsmi)
                    or checker.check_reaction("Reduction of ester to primary alcohol", rsmi)
                    or checker.check_reaction("Reduction of ketone to secondary alcohol", rsmi)
                    or checker.check_reaction(
                        "Reduction of carboxylic acid to primary alcohol", rsmi
                    )
                    or checker.check_reaction("Reduction of nitrile to amine", rsmi)
                    or checker.check_reaction("Reduction of primary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of secondary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of tertiary amides to amines", rsmi)
                    or checker.check_reaction("Reduction of nitro groups to amines", rsmi)
                ):
                    reduction_steps.append(depth)
                    print(f"Detected reduction at depth {depth}: {rsmi}")

                # If no specific reaction type matched, check for functional group changes
                if not (depth in oxidation_steps or depth in reduction_steps):
                    # Check for alcohol oxidation by looking at functional group changes
                    alcohol_in_product = (
                        checker.check_fg("Primary alcohol", product_smiles)
                        or checker.check_fg("Secondary alcohol", product_smiles)
                        or checker.check_fg("Tertiary alcohol", product_smiles)
                    )

                    aldehyde_in_reactants = any(
                        checker.check_fg("Aldehyde", r) for r in reactants_smiles
                    )
                    ketone_in_reactants = any(
                        checker.check_fg("Ketone", r) for r in reactants_smiles
                    )
                    carboxylic_in_reactants = any(
                        checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
                    )

                    if alcohol_in_product and (
                        aldehyde_in_reactants or ketone_in_reactants or carboxylic_in_reactants
                    ):
                        reduction_steps.append(depth)
                        print(f"Detected reduction (FG change) at depth {depth}: {rsmi}")

                    # Check for alcohol to carbonyl oxidation
                    alcohol_in_reactants = any(
                        checker.check_fg("Primary alcohol", r)
                        or checker.check_fg("Secondary alcohol", r)
                        or checker.check_fg("Tertiary alcohol", r)
                        for r in reactants_smiles
                    )

                    aldehyde_in_product = checker.check_fg("Aldehyde", product_smiles)
                    ketone_in_product = checker.check_fg("Ketone", product_smiles)
                    carboxylic_in_product = checker.check_fg("Carboxylic acid", product_smiles)

                    if alcohol_in_reactants and (
                        aldehyde_in_product or ketone_in_product or carboxylic_in_product
                    ):
                        oxidation_steps.append(depth)
                        print(f"Detected oxidation (FG change) at depth {depth}: {rsmi}")
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Oxidation steps found at depths: {oxidation_steps}")
    print(f"Reduction steps found at depths: {reduction_steps}")

    # In retrosynthesis, we're looking for:
    # Late-stage oxidation (lower depth) → middle-stage reduction → early-stage oxidation (higher depth)

    # Check for redox dance pattern in retrosynthetic context
    # We need to find: oxidation (low depth) → reduction (higher depth) → oxidation (even higher depth)
    for ox1 in oxidation_steps:
        for red in reduction_steps:
            if red > ox1:  # Reduction occurs earlier in synthesis (higher depth in retrosynthesis)
                for ox2 in oxidation_steps:
                    if ox2 > red:  # Second oxidation occurs even earlier (even higher depth)
                        print(
                            f"Redox dance pattern detected: oxidation at depth {ox1} → reduction at depth {red} → oxidation at depth {ox2}"
                        )
                        return True

    # Also check for the pattern in the opposite direction
    # This handles cases where the synthesis might be presented in forward direction
    for ox1 in oxidation_steps:
        for red in reduction_steps:
            if red < ox1:  # Reduction occurs later in synthesis (lower depth)
                for ox2 in oxidation_steps:
                    if ox2 < red:  # Second oxidation occurs even later (even lower depth)
                        print(
                            f"Redox dance pattern detected (forward direction): oxidation at depth {ox1} → reduction at depth {red} → oxidation at depth {ox2}"
                        )
                        return True

    # Check if we have at least one oxidation and one reduction step
    if oxidation_steps and reduction_steps:
        # Check if we have a pattern where oxidation and reduction alternate
        steps = [(depth, "oxidation") for depth in oxidation_steps] + [
            (depth, "reduction") for depth in reduction_steps
        ]
        steps.sort()  # Sort by depth

        # Check for alternating pattern
        for i in range(len(steps) - 2):
            if (
                steps[i][1] == "oxidation"
                and steps[i + 1][1] == "reduction"
                and steps[i + 2][1] == "oxidation"
            ) or (
                steps[i][1] == "reduction"
                and steps[i + 1][1] == "oxidation"
                and steps[i + 2][1] == "reduction"
            ):
                print(
                    f"Redox dance pattern detected (alternating): {steps[i]} → {steps[i+1]} → {steps[i+2]}"
                )
                return True

    return False
