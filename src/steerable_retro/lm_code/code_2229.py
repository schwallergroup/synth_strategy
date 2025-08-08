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
    Detects if the synthesis route includes a late-stage reductive amination.
    Late stage is defined as occurring in the first half of the synthesis depth.
    """
    max_depth = 0
    reductive_amination_depths = []

    # First pass to determine maximum depth
    def get_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            get_max_depth(child, current_depth + 1)

    # Second pass to detect reductive amination and its depth
    def detect_reductive_amination(node, current_depth=0):
        nonlocal reductive_amination_depths

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                print(f"Checking reaction at depth {current_depth}: {rsmi}")

                # Check for reductive amination using proper reaction checkers
                if (
                    checker.check_reaction("Reductive amination with aldehyde", rsmi)
                    or checker.check_reaction("Reductive amination with ketone", rsmi)
                    or checker.check_reaction("Reductive amination with alcohol", rsmi)
                ):
                    reductive_amination_depths.append(current_depth)
                    print(
                        f"Detected reductive amination at depth {current_depth} using reaction checker"
                    )
                else:
                    # Additional pattern recognition for reductive amination
                    try:
                        reactants_part = rsmi.split(">")[0]
                        products_part = rsmi.split(">")[-1]

                        reactants = reactants_part.split(".")
                        product = products_part

                        # Check if any reactant has a carbonyl group and another has an amine
                        has_carbonyl = any(
                            checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                            for r in reactants
                        )
                        has_amine = any(
                            checker.check_fg("Primary amine", r)
                            or checker.check_fg("Secondary amine", r)
                            for r in reactants
                        )

                        # Check if product has a new amine bond (where carbonyl was)
                        product_has_amine = (
                            checker.check_fg("Primary amine", product)
                            or checker.check_fg("Secondary amine", product)
                            or checker.check_fg("Tertiary amine", product)
                        )

                        # Look for specific pattern: nitrile reduction to amine
                        has_nitrile = any(checker.check_fg("Nitrile", r) for r in reactants)
                        has_carbonyl_compound = any(
                            checker.check_fg("Aldehyde", r) or checker.check_fg("Ketone", r)
                            for r in reactants
                        )

                        if (has_nitrile and has_carbonyl_compound and product_has_amine) or (
                            has_carbonyl and has_amine and product_has_amine
                        ):
                            reductive_amination_depths.append(current_depth)
                            print(
                                f"Detected reductive amination at depth {current_depth} using pattern recognition"
                            )
                    except Exception as e:
                        print(f"Error analyzing reaction: {e}")

        for child in node.get("children", []):
            detect_reductive_amination(child, current_depth + 1)

    # Run both passes
    get_max_depth(route)
    detect_reductive_amination(route)

    # Check if any reductive amination is in the first half of the synthesis (late stage)
    if reductive_amination_depths:
        # Find the earliest (lowest depth) reductive amination
        earliest_reductive_amination = min(reductive_amination_depths)
        # In retrosynthesis, lower depth means later stage
        is_late_stage = earliest_reductive_amination <= (max_depth // 2)
        print(
            f"Reductive amination is {'late' if is_late_stage else 'early'} stage (depth {earliest_reductive_amination}/{max_depth})"
        )
        return is_late_stage

    print("No reductive amination found in the route")
    return False
