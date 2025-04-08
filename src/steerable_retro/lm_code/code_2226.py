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
    This function detects if a synthesis preserves a heterocyclic core structure
    (specifically benzo[b]thiophene) throughout the synthesis.
    """
    # Track if the benzothiophene core is preserved throughout the synthesis
    core_preserved = True

    # Check if the final product (depth 0) has the benzothiophene core
    if route["type"] == "mol" and "smiles" in route:
        final_product_has_core = checker.check_ring("benzothiophene", route["smiles"])
        if not final_product_has_core:
            print("Final product does not contain benzothiophene core")
            return False

    # Function to check if the core is preserved in a reaction
    def check_reaction_preserves_core(reaction_node):
        if "metadata" not in reaction_node or "rsmi" not in reaction_node["metadata"]:
            return True  # Skip if no reaction SMILES available

        rsmi = reaction_node["metadata"]["rsmi"]
        reactants = rsmi.split(">")[0].split(".")
        product = rsmi.split(">")[-1]

        # Check if product has the core
        product_has_core = checker.check_ring("benzothiophene", product)

        # Check if any reactant has the core
        reactant_has_core = any(checker.check_ring("benzothiophene", r) for r in reactants)

        # Core is preserved if it's in both product and at least one reactant
        return product_has_core == reactant_has_core

    # Traverse the synthesis route
    def dfs_traverse(node):
        nonlocal core_preserved

        # If this is a reaction node, check if it preserves the core
        if node["type"] == "reaction" and not check_reaction_preserves_core(node):
            print(
                f"Core not preserved in reaction: {node.get('metadata', {}).get('rsmi', 'unknown')}"
            )
            core_preserved = False
            return

        # Process children
        for child in node.get("children", []):
            if core_preserved:  # Only continue if core is still preserved
                dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    if core_preserved:
        print("Heterocyclic core (benzothiophene) is preserved throughout synthesis")
    else:
        print("Heterocyclic core is not preserved throughout synthesis")

    return core_preserved
