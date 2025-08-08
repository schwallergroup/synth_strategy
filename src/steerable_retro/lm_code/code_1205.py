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
    Detects if the synthetic route involves late-stage sulfonamide formation.
    Looks for sulfonamide formation in the second half of the synthesis.
    """
    max_depth = 0
    sulfonamide_depth = None

    # First pass to find max depth
    def find_max_depth(node, current_depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, current_depth)

        for child in node.get("children", []):
            find_max_depth(child, current_depth + 1)

    # Second pass to find sulfonamide formation
    def find_sulfonamide(node, depth=0):
        nonlocal sulfonamide_depth

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            print(f"Checking reaction at depth {depth}: {rsmi}")

            # Check if reactants contain sulfonyl halide
            reactants = reactants_str.split(".")
            has_sulfonyl_halide = any(checker.check_fg("Sulfonyl halide", r) for r in reactants)
            print(f"  Has sulfonyl halide in reactants: {has_sulfonyl_halide}")

            # Check if product has sulfonamide pattern
            has_sulfonamide = checker.check_fg("Sulfonamide", product_str)
            print(f"  Has sulfonamide in product: {has_sulfonamide}")

            # Check if reactants do NOT have sulfonamide (to ensure it's being formed)
            no_sulfonamide_in_reactants = not any(
                checker.check_fg("Sulfonamide", r) for r in reactants
            )
            print(f"  No sulfonamide in reactants: {no_sulfonamide_in_reactants}")

            # Check if this is a sulfonamide formation reaction
            is_sulfonamide_reaction = (
                checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                )
                or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )
                or checker.check_reaction("Schotten-Baumann to ester", rsmi)  # Sometimes mislabeled
            )
            print(f"  Is recognized sulfonamide reaction: {is_sulfonamide_reaction}")

            # Alternative check: look for sulfonyl halide + amine reaction pattern
            has_amine = any(
                checker.check_fg("Primary amine", r)
                or checker.check_fg("Secondary amine", r)
                or checker.check_fg("Aniline", r)
                for r in reactants
            )
            print(f"  Has amine in reactants: {has_amine}")

            # Determine if this is a sulfonamide formation
            is_forming_sulfonamide = (
                has_sulfonamide
                and has_sulfonyl_halide
                and has_amine
                and no_sulfonamide_in_reactants
                and (is_sulfonamide_reaction or (has_sulfonyl_halide and has_amine))
            )

            if is_forming_sulfonamide:
                print(f"Detected sulfonamide formation at depth {depth}")
                # Track the earliest (lowest depth) sulfonamide formation
                if sulfonamide_depth is None or depth < sulfonamide_depth:
                    sulfonamide_depth = depth

        for child in node.get("children", []):
            find_sulfonamide(child, depth + 1)

    find_max_depth(route)
    print(f"Maximum synthesis depth: {max_depth}")
    find_sulfonamide(route)

    # Check if sulfonamide formation occurs in second half of synthesis
    if sulfonamide_depth is not None:
        is_late_stage = sulfonamide_depth <= (max_depth // 2)
        print(f"Sulfonamide formation at depth {sulfonamide_depth} of max depth {max_depth}")
        print(f"Is late stage: {is_late_stage}")
        return is_late_stage
    else:
        print("No sulfonamide formation detected in the route")

    return False
