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
    This function detects a synthesis strategy involving isoxazole heterocycle disconnection.
    In retrosynthesis, we're looking for reactions where isoxazole is formed in the forward direction,
    which means it should be in the products but not in the reactants.
    """
    isoxazole_disconnection_found = False

    def dfs_traverse(node):
        nonlocal isoxazole_disconnection_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            # Extract reaction SMILES
            rsmi = node["metadata"]["rsmi"]

            # Split into reactants and product
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if any reactant contains isoxazole
            isoxazole_in_reactants = False
            for reactant in reactants:
                if checker.check_ring("isoxazole", reactant):
                    isoxazole_in_reactants = True
                    break

            # Check if product contains isoxazole
            isoxazole_in_products = checker.check_ring("isoxazole", product)

            # Check for isoxazole disconnection: present in products but not in reactants
            # This means isoxazole is formed in the forward reaction
            if isoxazole_in_products and not isoxazole_in_reactants:
                print(
                    f"Found isoxazole disconnection: isoxazole in product but not in reactants. Reaction: {rsmi}"
                )
                isoxazole_disconnection_found = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return isoxazole_disconnection_found
