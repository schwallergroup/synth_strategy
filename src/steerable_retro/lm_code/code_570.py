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
    This function detects if the final step in the synthesis is a Boc deprotection.
    In a retrosynthetic tree, the final synthetic step would be at the leaf reaction nodes.
    """
    # Define Boc deprotection reaction types
    boc_deprotection_types = [
        "Boc amine deprotection",
        "Boc amine deprotection of guanidine",
        "Boc amine deprotection to NH-NH2",
        "Tert-butyl deprotection of amine",
    ]

    # Helper function to check if a reaction is a Boc deprotection
    def is_boc_deprotection(reaction_node):
        if not reaction_node.get("metadata", {}).get("rsmi"):
            return False

        rsmi = reaction_node["metadata"]["rsmi"]

        # Check if the reaction is a Boc deprotection using the checker function
        for boc_type in boc_deprotection_types:
            if checker.check_reaction(boc_type, rsmi):
                print(f"Found {boc_type} as final step: {rsmi}")
                return True

        # If none of the specific Boc deprotection reactions match,
        # check manually for Boc group removal
        reactants = rsmi.split(">")[0].split(".")
        product = rsmi.split(">")[-1]

        # Check for Boc group in reactants but not in product
        reactant_has_boc = any(checker.check_fg("Boc", r) for r in reactants if r)
        product_has_boc = checker.check_fg("Boc", product) if product else False

        if reactant_has_boc and not product_has_boc:
            print(f"Found Boc deprotection as final step (manual check): {rsmi}")
            return True

        return False

    # Find all leaf reaction nodes (final steps in forward synthesis)
    found_boc_deprotection = [False]  # Using list to allow modification in nested function

    def find_leaf_reactions(node, is_child_of_reaction=False):
        if node["type"] == "reaction":
            # Check if this reaction node has only molecule children (leaf reaction)
            all_children_are_mols = all(
                child["type"] == "mol" for child in node.get("children", [])
            )

            if all_children_are_mols:
                # This is a leaf reaction node (final step in forward synthesis)
                if is_boc_deprotection(node):
                    found_boc_deprotection[0] = True
                    return

            # Continue traversal for non-leaf reaction nodes
            for child in node.get("children", []):
                find_leaf_reactions(child, True)

        elif node["type"] == "mol" and not is_child_of_reaction:
            # For molecule nodes that aren't children of reaction nodes
            for child in node.get("children", []):
                find_leaf_reactions(child, False)

    # Start traversal from the root
    find_leaf_reactions(route)

    return found_boc_deprotection[0]
