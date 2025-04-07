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
    This function detects routes with both pyrazole and oxazole heterocycles in the final product.
    """
    # In a retrosynthetic tree, the final product is the root node
    if route["type"] != "mol":
        print("Root node is not a molecule, cannot check for heterocycles")
        return False

    # Check if the final product contains both pyrazole and oxazole
    smiles = route["smiles"]
    print(f"Checking final product: {smiles}")

    has_pyrazole = checker.check_ring("pyrazole", smiles)
    has_oxazole = checker.check_ring("oxazole", smiles)

    print(f"Has pyrazole: {has_pyrazole}")
    print(f"Has oxazole: {has_oxazole}")

    if has_pyrazole and has_oxazole:
        print(f"Found both pyrazole and oxazole in final product: {smiles}")
        return True
    else:
        print(f"Did not find both pyrazole and oxazole in the final product")
        return False
