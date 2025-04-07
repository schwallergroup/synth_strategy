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
    This function detects if the synthetic route involves a late-stage condensation to form an oxime.
    Oximes are typically formed by the reaction of a carbonyl compound (aldehyde or ketone) with hydroxylamine.
    """
    late_stage_oxime_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal late_stage_oxime_formation

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            # Check if product contains oxime
            product_has_oxime = checker.check_fg("Oxime", product_part)

            if product_has_oxime and depth <= 1:  # Late stage = depth 0 or 1
                print(f"Found product with oxime at depth {depth}: {product_part}")

                # Check if reactants don't contain oxime (meaning it was formed in this reaction)
                reactants = reactants_part.split(".")
                reactants_have_oxime = any(
                    checker.check_fg("Oxime", reactant) for reactant in reactants
                )

                if not reactants_have_oxime:
                    print(f"Oxime not present in reactants - potential formation reaction")

                    # Check for carbonyl groups in reactants (needed for oxime formation)
                    has_carbonyl = any(
                        checker.check_fg("Aldehyde", reactant)
                        or checker.check_fg("Ketone", reactant)
                        or checker.check_fg("Formaldehyde", reactant)
                        for reactant in reactants
                    )

                    # Improved hydroxylamine detection
                    has_hydroxylamine = any(
                        "NOH" in reactant
                        or "NH2OH" in reactant
                        or "N(O)" in reactant
                        or "NO" in reactant
                        for reactant in reactants
                    )

                    # Check if this is a condensation reaction - expanded checks
                    is_condensation = (
                        checker.check_reaction(
                            "Addition of primary amines to aldehydes/thiocarbonyls", rsmi
                        )
                        or checker.check_reaction(
                            "Addition of primary amines to ketones/thiocarbonyls", rsmi
                        )
                        or checker.check_reaction("Ketone/aldehyde to hydrazone", rsmi)
                        or checker.check_reaction("Aldol condensation", rsmi)
                    )

                    # Pattern-based detection for oxime formation
                    # Look for C=O → C=N-OH transformation pattern
                    oxime_formation_pattern = False
                    if has_carbonyl and not is_condensation:
                        # If we have carbonyl in reactants and oxime in product, it's likely an oxime formation
                        # even if not explicitly matched by reaction patterns
                        oxime_formation_pattern = True

                    print(
                        f"Has carbonyl: {has_carbonyl}, Has hydroxylamine: {has_hydroxylamine}, Is condensation: {is_condensation}, Oxime pattern: {oxime_formation_pattern}"
                    )

                    # If we have the right conditions for oxime formation at a late stage
                    if has_carbonyl and (
                        has_hydroxylamine or is_condensation or oxime_formation_pattern
                    ):
                        print(f"Late-stage oxime formation detected in reaction: {rsmi}")
                        late_stage_oxime_formation = True

        # Traverse children with incremented depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    print(f"Late-stage condensation strategy detected: {late_stage_oxime_formation}")
    return late_stage_oxime_formation
