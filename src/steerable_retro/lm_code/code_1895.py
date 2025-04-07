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
    This function detects a late-stage N-dealkylation, specifically the removal
    of an allyl group from a tertiary amine.
    """

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Only check reactions at depth 0 or 1 (late stage)
            if depth <= 1:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                print(f"Checking reaction at depth {depth}: {rsmi}")

                # In forward reaction: Tertiary amine -> Secondary amine (N-dealkylation)
                # In retrosynthesis: Secondary amine -> Tertiary amine

                # Check if reactant has tertiary amine (in forward direction)
                reactant_has_tertiary_amine = False
                tertiary_amine_reactant = None
                for r in reactants:
                    if checker.check_fg("Tertiary amine", r):
                        reactant_has_tertiary_amine = True
                        tertiary_amine_reactant = r
                        break

                if not reactant_has_tertiary_amine:
                    print("No reactant with tertiary amine found")
                    return False

                # Check if product has secondary amine
                if not checker.check_fg("Secondary amine", product):
                    print("Product does not have secondary amine")
                    return False

                # Check if reactant has allyl group
                if not checker.check_fg("Allyl", tertiary_amine_reactant):
                    print("Tertiary amine reactant does not have allyl group")
                    return False

                # Check if this is a dealkylation reaction
                if checker.check_reaction("Hydrogenolysis of tertiary amines", rsmi):
                    print(f"Found late-stage N-dealkylation via hydrogenolysis at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    print(f"Reactant (tertiary amine with allyl): {tertiary_amine_reactant}")
                    print(f"Product (secondary amine): {product}")
                    return True

                # Alternative check for N-dealkylation
                # Look for a reaction where a tertiary amine loses an allyl group
                # This is a more general check for N-dealkylation
                print("Checking for general N-dealkylation pattern")

                # If we have a tertiary amine in reactant, secondary amine in product,
                # and allyl in reactant but not in product, it's likely N-dealkylation
                if not checker.check_fg("Allyl", product) or (
                    checker.check_fg("Allyl", product)
                    and checker.check_fg("Allyl", tertiary_amine_reactant)
                    and tertiary_amine_reactant.count("C=CC") > product.count("C=CC")
                ):

                    print(f"Found late-stage N-dealkylation (general pattern) at depth {depth}")
                    print(f"Reaction SMILES: {rsmi}")
                    print(f"Reactant (tertiary amine with allyl): {tertiary_amine_reactant}")
                    print(f"Product (secondary amine): {product}")
                    return True

        # Recursively check children
        for child in node.get("children", []):
            if dfs_traverse(child, depth + 1):
                return True

        return False

    return dfs_traverse(route)
