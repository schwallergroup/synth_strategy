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
    This function detects if a halogen substituent is maintained throughout the synthesis.
    """
    # Track halogen atoms through the synthesis
    halogen_maintained = False

    # Get the final product (root node)
    final_product_smiles = route["smiles"]
    final_product = Chem.MolFromSmiles(final_product_smiles)

    # Check if final product has any halogen
    has_halogen = False
    halogen_types = [
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Aromatic halide",
        "Alkenyl halide",
        "Haloalkyne",
    ]

    for halogen_type in halogen_types:
        if checker.check_fg(halogen_type, final_product_smiles):
            has_halogen = True
            print(f"Final product contains {halogen_type}")
            break

    if not has_halogen:
        print("Final product does not contain any halogen")
        return False

    # Track if halogen is maintained throughout synthesis
    def dfs_traverse(node, depth=0):
        nonlocal halogen_maintained

        # Process molecule nodes
        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Skip checking starting materials (in_stock)
            if node.get("in_stock", False):
                return

            # Check if this intermediate has a halogen
            has_halogen_intermediate = False
            for halogen_type in halogen_types:
                if checker.check_fg(halogen_type, mol_smiles):
                    has_halogen_intermediate = True
                    print(f"Intermediate at depth {depth} contains {halogen_type}")
                    break

            if not has_halogen_intermediate and depth > 0:  # Skip final product (depth 0)
                print(f"Halogen not maintained in intermediate at depth {depth}")
                halogen_maintained = False
                return

        # Process reaction nodes and check if halogen is preserved
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if this reaction preserves halogens
            # Skip reactions that are known to remove or modify halogens
            halogen_removing_reactions = ["Aromatic dehalogenation", "Dehalogenation"]

            for rxn_type in halogen_removing_reactions:
                if checker.check_reaction(rxn_type, rsmi):
                    print(f"Halogen-removing reaction detected: {rxn_type}")
                    halogen_maintained = False
                    return

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal and set initial value
    halogen_maintained = True
    dfs_traverse(route)

    if halogen_maintained:
        print("Halogen substituent maintained throughout synthesis")
    else:
        print("Halogen substituent NOT maintained throughout synthesis")

    return halogen_maintained
