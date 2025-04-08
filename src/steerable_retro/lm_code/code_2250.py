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
    This function detects if the synthesis involves a late-stage
    nitro reduction to form an amine.
    """
    # Track if nitro reduction occurs at late stage (depth ≤ 2)
    late_stage_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_reduction

        if node["type"] == "reaction" and depth <= 2:  # Late stage (depth 0-2)
            # Extract reactants and product
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    return

                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Primary check: Use the reaction checker
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print(f"Found nitro reduction reaction at depth {depth}: {rsmi}")
                    late_stage_reduction = True
                    return

                # Alternative check: verify nitro group in reactant and amine in product
                product_has_primary_amine = checker.check_fg("Primary amine", product_smiles)
                product_has_secondary_amine = checker.check_fg("Secondary amine", product_smiles)
                product_has_tertiary_amine = checker.check_fg("Tertiary amine", product_smiles)
                product_has_aniline = checker.check_fg("Aniline", product_smiles)

                product_has_any_amine = (
                    product_has_primary_amine
                    or product_has_secondary_amine
                    or product_has_tertiary_amine
                    or product_has_aniline
                )

                if product_has_any_amine:
                    # Check if any reactant has nitro group
                    for reactant in reactants_smiles:
                        if checker.check_fg("Nitro group", reactant):
                            print(
                                f"Found nitro group in reactant and amine in product at depth {depth}"
                            )

                            # Additional verification: check if product doesn't have nitro group
                            if not checker.check_fg("Nitro group", product_smiles):
                                print(
                                    f"Confirmed: Nitro group was reduced to amine at depth {depth}"
                                )
                                late_stage_reduction = True
                                return
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late stage nitro reduction detected: {late_stage_reduction}")
    return late_stage_reduction
