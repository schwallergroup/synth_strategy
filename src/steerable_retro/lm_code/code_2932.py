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
    This function detects if the synthesis uses late-stage sulfonamide formation
    as the final step in the synthesis.
    """
    final_step_is_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_is_sulfonamide

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction" and depth <= 1:  # Final or near-final step
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check if this is a sulfonamide formation reaction
                is_primary_sulfonamide = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rsmi
                )
                is_secondary_sulfonamide = checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rsmi
                )

                # Check if product contains sulfonamide
                product_has_sulfonamide = checker.check_fg("Sulfonamide", product)
                print(f"Product has sulfonamide: {product_has_sulfonamide}")

                # Check if any reactant doesn't have sulfonamide
                reactant_without_sulfonamide = False
                for reactant in reactants:
                    if not checker.check_fg("Sulfonamide", reactant):
                        reactant_without_sulfonamide = True
                        break

                print(f"Reactant without sulfonamide: {reactant_without_sulfonamide}")
                print(f"Is primary sulfonamide reaction: {is_primary_sulfonamide}")
                print(f"Is secondary sulfonamide reaction: {is_secondary_sulfonamide}")

                # Check for required functional groups in reactants
                sulfonyl_chloride_present = any(
                    checker.check_fg("Sulfonyl halide", r) for r in reactants
                )
                amine_present = any(
                    checker.check_fg("Primary amine", r) or checker.check_fg("Secondary amine", r)
                    for r in reactants
                )

                print(f"Sulfonyl chloride present: {sulfonyl_chloride_present}")
                print(f"Amine present: {amine_present}")

                # Check for sulfonamide formation - either by named reaction or by functional group presence
                if product_has_sulfonamide and reactant_without_sulfonamide and depth <= 1:
                    if is_primary_sulfonamide or is_secondary_sulfonamide:
                        print(f"Found late-stage sulfonamide formation (named reaction): {rsmi}")
                        final_step_is_sulfonamide = True
                    elif sulfonyl_chloride_present and amine_present:
                        print(f"Found late-stage sulfonamide formation (functional groups): {rsmi}")
                        final_step_is_sulfonamide = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return final_step_is_sulfonamide
