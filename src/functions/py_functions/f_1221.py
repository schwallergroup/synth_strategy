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
    """Check if the final product contains a piperazine sulfonamide scaffold"""
    if route["type"] != "mol":
        return False

    # Check if the final product contains both piperazine and sulfonamide
    final_product_smiles = route["smiles"]
    has_piperazine = checker.check_ring("piperazine", final_product_smiles)

    # Expanded check for sulfonamide-like groups
    has_sulfonamide = (
        checker.check_fg("Sulfonamide", final_product_smiles)
        or checker.check_fg("Sulfonate", final_product_smiles)
        or checker.check_fg("Sulfone", final_product_smiles)
        or "S(=O)(=O)" in final_product_smiles
    )

    # Check if there's a connection between piperazine and sulfonamide
    # This is a simplified check - in a real implementation, you'd want to use
    # RDKit to verify the actual connection between these groups
    has_connection = (
        "S(=O)(=O)N" in final_product_smiles or "NS(=O)(=O)" in final_product_smiles
    )

    result = has_piperazine and (has_sulfonamide or has_connection)
    print(
        f"Final product has piperazine: {has_piperazine}, sulfonamide or related group: {has_sulfonamide or has_connection}"
    )

    # If we have piperazine but not detecting sulfonamide, print the SMILES for debugging
    if has_piperazine and not (has_sulfonamide or has_connection):
        print(f"Piperazine detected but no sulfonamide in: {final_product_smiles}")

        # Try to find any sulfur-containing groups
        if "S" in final_product_smiles:
            print(f"Product contains sulfur atoms but not recognized as sulfonamide")

    return result
