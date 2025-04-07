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
    Detects preservation of halogen (F, Cl, Br, I) throughout the synthesis.
    Checks if halogens are present in both starting materials and final product.
    """
    halogen_groups = [
        "Aromatic halide",
        "Primary halide",
        "Secondary halide",
        "Tertiary halide",
        "Alkenyl halide",
        "Haloalkyne",
    ]

    final_product_has_halogen = False
    starting_material_has_halogen = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_halogen, starting_material_has_halogen

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check if molecule contains any halogen
            has_halogen = any(checker.check_fg(fg, mol_smiles) for fg in halogen_groups)

            if has_halogen:
                if depth == 0:
                    # This is the final product (root node)
                    final_product_has_halogen = True
                    print(f"Found halogen in final product: {mol_smiles}")
                elif "in_stock" in node and node["in_stock"]:
                    # This is a starting material
                    starting_material_has_halogen = True
                    print(f"Found halogen in starting material: {mol_smiles}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return final_product_has_halogen and starting_material_has_halogen
