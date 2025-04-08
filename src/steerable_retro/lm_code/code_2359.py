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
    Detects late-stage sulfonylation strategy in the synthesis route.
    """
    sulfonylation_found = False

    def dfs(node, depth=0):
        nonlocal sulfonylation_found

        if node["type"] == "reaction":
            try:
                rxn_smiles = node.get("metadata", {}).get("rsmi", "")
                if not rxn_smiles:
                    return

                # Check for sulfonylation reactions
                if (
                    checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) primary amine", rxn_smiles
                    )
                    or checker.check_reaction(
                        "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rxn_smiles
                    )
                    or checker.check_reaction("Formation of Sulfonic Esters", rxn_smiles)
                    or checker.check_reaction(
                        "Formation of Sulfonic Esters on TMS protected alcohol", rxn_smiles
                    )
                    or checker.check_reaction("sulfon_amide", rxn_smiles)
                ):

                    # Check if this is late-stage (depth <= 3)
                    if depth <= 3:
                        # Verify that a sulfonyl group is actually formed
                        reactants = rxn_smiles.split(">")[0].split(".")
                        product = rxn_smiles.split(">")[-1]

                        # Check if product has a sulfonamide, sulfonate, or sulfone group
                        has_sulfonyl = (
                            checker.check_fg("Sulfonamide", product)
                            or checker.check_fg("Sulfonate", product)
                            or checker.check_fg("Sulfone", product)
                            or checker.check_fg("Sulfamate", product)
                            or checker.check_fg("Sulfonyl halide", product)
                        )

                        # Check if reactants don't already have these groups
                        reactants_have_sulfonyl = any(
                            checker.check_fg("Sulfonamide", r)
                            or checker.check_fg("Sulfonate", r)
                            or checker.check_fg("Sulfone", r)
                            or checker.check_fg("Sulfamate", r)
                            or checker.check_fg("Sulfonyl halide", r)
                            for r in reactants
                        )

                        if has_sulfonyl and not reactants_have_sulfonyl:
                            sulfonylation_found = True
                            print(
                                f"Found late-stage sulfonylation reaction at depth {depth}: {rxn_smiles}"
                            )
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    print(f"Late-stage sulfonylation found: {sulfonylation_found}")
    return sulfonylation_found
