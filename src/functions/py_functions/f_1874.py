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
    Detects a synthetic strategy based on an indole core scaffold that is maintained
    throughout the synthesis.
    """
    # Track if the final product contains indole or quinoline
    final_product_has_indole = False

    # Track reactions where indole/quinoline core is maintained
    reactions_maintaining_indole = 0
    total_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_indole, reactions_maintaining_indole, total_reactions

        # Check if current node is a molecule
        if node["type"] == "mol" and depth == 0:
            # This is the final product (depth 0)
            mol_smiles = node["smiles"]
            # Check for indole or quinoline structure
            final_product_has_indole = (
                checker.check_ring("indole", mol_smiles)
                or checker.check_ring("quinoline", mol_smiles)
                or checker.check_ring("isoquinoline", mol_smiles)
            )
            print(
                f"Final product has indole/quinoline: {final_product_has_indole}, SMILES: {mol_smiles}"
            )

        # Check if current node is a reaction
        elif node["type"] == "reaction":
            total_reactions += 1
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if product contains indole or quinoline structure
            product_has_indole = (
                checker.check_ring("indole", product)
                or checker.check_ring("quinoline", product)
                or checker.check_ring("isoquinoline", product)
            )
            print(
                f"Product has indole/quinoline: {product_has_indole}, SMILES: {product}"
            )

            # Check if any reactant contains indole or quinoline structure
            reactants_have_indole = any(
                checker.check_ring("indole", reactant)
                or checker.check_ring("quinoline", reactant)
                or checker.check_ring("isoquinoline", reactant)
                for reactant in reactants
            )
            print(
                f"Reactants have indole/quinoline: {reactants_have_indole}, SMILES: {reactants}"
            )

            # If both reactants and product have indole/quinoline, the core is maintained
            if product_has_indole and reactants_have_indole:
                reactions_maintaining_indole += 1
                print(f"Indole/quinoline core maintained in reaction {total_reactions}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(
        f"Final stats: reactions_maintaining_indole={reactions_maintaining_indole}, total_reactions={total_reactions}"
    )

    # Return True if:
    # 1. The final product has indole/quinoline
    # 2. At least 3 reactions maintain the indole/quinoline core
    # 3. We have at least 3 total reactions
    return (
        final_product_has_indole
        and reactions_maintaining_indole >= 3
        and total_reactions >= 3
    )
