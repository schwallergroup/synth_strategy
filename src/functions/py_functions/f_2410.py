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
    Detects if a reaction in the synthetic route is an amide formation by checking for:
    1. Presence of carboxylic acid in reactants
    2. Presence of amine in reactants
    3. Formation of amide bond in product
    4. Reaction type matches amide formation patterns
    """

    def check_node(node):
        # Only process reaction nodes
        if node.get("type") != "reaction":
            return False

        try:
            # Extract reaction SMILES from metadata
            if "rsmi" not in node.get("metadata", {}):
                return False

            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants_smiles = reactants_part.split(".")
            product_smiles = product_part

            # Check if this is an amide formation reaction using the checker
            if checker.check_reaction(
                "Carboxylic acid with primary amine to amide", rsmi
            ):
                print(f"Found amide formation reaction: {rsmi}")
                return True

            if checker.check_reaction(
                "Acyl chloride with primary amine to amide (Schotten-Baumann)", rsmi
            ):
                print(f"Found Schotten-Baumann amide formation: {rsmi}")
                return True

            if checker.check_reaction("Ester with primary amine to amide", rsmi):
                print(f"Found ester-amine amide formation: {rsmi}")
                return True

            # If specific reaction checks failed, try a more general approach
            # Check for carboxylic acid in reactants
            has_carboxylic_acid = any(
                checker.check_fg("Carboxylic acid", r) for r in reactants_smiles
            )

            # Check for acyl chloride in reactants (alternative to carboxylic acid)
            has_acyl_chloride = any(
                checker.check_fg("Acyl halide", r) for r in reactants_smiles
            )

            # Check for ester in reactants (alternative to carboxylic acid)
            has_ester = any(checker.check_fg("Ester", r) for r in reactants_smiles)

            # Check for amine in reactants
            has_primary_amine = any(
                checker.check_fg("Primary amine", r) for r in reactants_smiles
            )
            has_secondary_amine = any(
                checker.check_fg("Secondary amine", r) for r in reactants_smiles
            )
            has_amine = has_primary_amine or has_secondary_amine

            # Check for amide in product
            has_amide_product = (
                checker.check_fg("Primary amide", product_smiles)
                or checker.check_fg("Secondary amide", product_smiles)
                or checker.check_fg("Tertiary amide", product_smiles)
            )

            # If reactants have carboxylic acid/acyl chloride/ester and amine, and product has amide, it's likely an amide formation
            if (
                (has_carboxylic_acid or has_acyl_chloride or has_ester)
                and has_amine
                and has_amide_product
            ):
                print(f"Detected amide formation based on functional groups: {rsmi}")
                print(
                    f"Carboxylic acid: {has_carboxylic_acid}, Acyl halide: {has_acyl_chloride}, Ester: {has_ester}"
                )
                print(f"Amine: {has_amine}, Amide in product: {has_amide_product}")
                return True

            return False

        except Exception as e:
            print(f"Error in amide formation detection: {e}")
            return False

    # Traverse the route to find reaction nodes
    def traverse(node):
        if check_node(node):
            return True

        # Check children nodes
        for child in node.get("children", []):
            if traverse(child):
                return True

        return False

    # Start traversal from the root
    return traverse(route)
