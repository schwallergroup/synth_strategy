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
    Detects if cyano and chloro functional groups are preserved throughout the synthesis.
    This indicates a strategy that builds around these sensitive functional groups.
    """
    final_product_has_cyano = False
    final_product_has_chloro = False
    starting_material_has_cyano = False
    starting_material_has_chloro = False

    def dfs_traverse(node, depth=0):
        nonlocal final_product_has_cyano, final_product_has_chloro
        nonlocal starting_material_has_cyano, starting_material_has_chloro

        if node["type"] == "mol":
            mol_smiles = node["smiles"]

            # Check for cyano and chloro groups using checker functions
            has_cyano = checker.check_fg("Nitrile", mol_smiles)
            has_chloro = (
                checker.check_fg("Primary halide", mol_smiles)
                or checker.check_fg("Secondary halide", mol_smiles)
                or checker.check_fg("Tertiary halide", mol_smiles)
                or checker.check_fg("Aromatic halide", mol_smiles)
            )

            if depth == 0:  # Final product
                final_product_has_cyano = has_cyano
                final_product_has_chloro = has_chloro
                if has_cyano:
                    print(f"Final product has cyano group: {mol_smiles}")
                if has_chloro:
                    print(f"Final product has chloro group: {mol_smiles}")

            if node.get("in_stock", False):  # Starting material
                if has_cyano:
                    starting_material_has_cyano = True
                    print(f"Starting material has cyano group: {mol_smiles}")
                if has_chloro:
                    starting_material_has_chloro = True
                    print(f"Starting material has chloro group: {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Strategy is present if both functional groups are preserved from starting materials to final product
    result = (
        final_product_has_cyano
        and starting_material_has_cyano
        and final_product_has_chloro
        and starting_material_has_chloro
    )

    print(f"Final result: {result}")
    print(
        f"Final product has cyano: {final_product_has_cyano}, Starting material has cyano: {starting_material_has_cyano}"
    )
    print(
        f"Final product has chloro: {final_product_has_chloro}, Starting material has chloro: {starting_material_has_chloro}"
    )

    return result
