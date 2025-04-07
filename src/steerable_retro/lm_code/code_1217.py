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
    """Check for sulfonamide formation in the synthesis route"""
    found = False

    def dfs(node, depth=0):
        nonlocal found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rxn_smiles = node["metadata"]["rsmi"]

            # Check for sulfonamide formation reactions
            if (
                checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) primary amine", rxn_smiles
                )
                or checker.check_reaction(
                    "Sulfonamide synthesis (Schotten-Baumann) secondary amine", rxn_smiles
                )
                or checker.check_reaction("sulfon_amide", rxn_smiles)
            ):

                # Verify that sulfonamide is formed
                try:
                    product = rxn_smiles.split(">")[-1]
                    if (
                        checker.check_fg("Sulfonamide", product)
                        or checker.check_fg("Sulfonate", product)
                        or checker.check_fg("Sulfone", product)
                    ):
                        found = True
                        print(f"Found sulfonamide formation at depth {depth}: {rxn_smiles}")
                except Exception as e:
                    print(f"Error checking sulfonamide formation: {e}")

            # Additional check for sulfonamide formation even if not explicitly labeled
            if not found:
                try:
                    reactants = rxn_smiles.split(">")[0].split(".")
                    product = rxn_smiles.split(">")[-1]

                    has_sulfonyl = any(checker.check_fg("Sulfonyl halide", r) for r in reactants)
                    has_amine = any(
                        checker.check_fg("Primary amine", r)
                        or checker.check_fg("Secondary amine", r)
                        or checker.check_fg("Aniline", r)
                        for r in reactants
                    )
                    has_sulfonamide = (
                        checker.check_fg("Sulfonamide", product)
                        or checker.check_fg("Sulfonate", product)
                        or checker.check_fg("Sulfone", product)
                    )

                    if has_sulfonyl and has_amine and has_sulfonamide:
                        found = True
                        print(
                            f"Found implicit sulfonamide formation at depth {depth}: {rxn_smiles}"
                        )
                except Exception as e:
                    print(f"Error checking implicit sulfonamide formation: {e}")

        # Recursively check children
        for child in node.get("children", []):
            dfs(child, depth + 1)

    dfs(route)
    return found
