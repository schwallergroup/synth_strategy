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
    This function detects if an indole scaffold is preserved throughout the synthesis route.
    """
    indole_present_in_all_products = True

    def dfs_traverse(node, depth=0):
        nonlocal indole_present_in_all_products

        if not indole_present_in_all_products:
            return  # Early termination if we already found a violation

        if node["type"] == "mol" and not node.get("in_stock", False):
            # Check if this molecule contains an indole scaffold
            mol_smiles = node["smiles"]
            if not checker.check_ring("indole", mol_smiles):
                indole_present_in_all_products = False
                print(f"Molecule without indole scaffold found: {mol_smiles}")
                return

        elif node["type"] == "reaction":
            # For reaction nodes, check if the indole scaffold is preserved
            try:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has indole
                if not checker.check_ring("indole", product_smiles):
                    indole_present_in_all_products = False
                    print(f"Product without indole scaffold found in reaction: {product_smiles}")
                    return

                # Check if at least one reactant has indole
                # This ensures the indole scaffold is not created in this reaction
                indole_in_reactants = any(checker.check_ring("indole", r) for r in reactants_smiles)
                if not indole_in_reactants:
                    indole_present_in_all_products = False
                    print(
                        f"No reactant with indole scaffold found in reaction. Indole was created, not preserved."
                    )
                    return

            except (KeyError, IndexError) as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    print(f"Indole scaffold preservation: {indole_present_in_all_products}")
    return indole_present_in_all_products
