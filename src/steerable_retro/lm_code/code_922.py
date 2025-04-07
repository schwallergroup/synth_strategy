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
    This function detects if the final step in a synthesis involves
    oxidation of a sulfide to a sulfoxide.
    """
    found_late_stage_oxidation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_late_stage_oxidation

        # The final product is at depth 0 (mol node)
        # The final reaction step is at depth 1 (reaction node)
        if (
            depth == 1
            and node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0]
            product_smiles = rsmi.split(">")[-1]

            try:
                # Check if this is a sulfide to sulfoxide oxidation
                if checker.check_reaction("Sulfanyl to sulfinyl", rsmi):
                    print(f"Found late-stage sulfide oxidation via reaction check: {rsmi}")
                    found_late_stage_oxidation = True
                else:
                    # Alternative check: verify sulfide in reactants becomes sulfoxide in product
                    reactants = reactants_smiles.split(".")

                    # Check if any reactant has a sulfide but not a sulfoxide
                    has_sulfide_reactant = any(
                        checker.check_fg("Monosulfide", r) and not checker.check_fg("Sulfoxide", r)
                        for r in reactants
                    )

                    # Check if product has a sulfoxide
                    has_sulfoxide_product = checker.check_fg("Sulfoxide", product_smiles)

                    if has_sulfide_reactant and has_sulfoxide_product:
                        print(f"Found late-stage sulfide oxidation via FG check: {rsmi}")
                        found_late_stage_oxidation = True
            except Exception as e:
                print(f"Error processing reaction at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Late-stage sulfoxide formation detection result: {found_late_stage_oxidation}")
    return found_late_stage_oxidation
