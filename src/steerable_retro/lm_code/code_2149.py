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
    This function detects if the final product contains multiple heterocycles
    (including indole, indazole, tetrazole, oxadiazole, and others).
    """
    # First, find the final product (root node)
    final_product = route["smiles"] if route["type"] == "mol" else None

    if not final_product:
        return False

    try:
        # Define heterocycles to check
        heterocycles = [
            "indole",
            "indazole",
            "tetrazole",
            "oxadiazole",
            "benzimidazole",
            "benzoxazole",
            "benzothiazole",
            "triazole",
            "pyrazole",
            "imidazole",
            "thiazole",
            "oxazole",
            "isoxazole",
            "pyridine",
            "pyrimidine",
            "pyrazine",
            "pyridazine",
        ]

        # Count the number of different heterocycles
        heterocycle_count = 0
        for heterocycle in heterocycles:
            if checker.check_ring(heterocycle, final_product):
                heterocycle_count += 1
                print(f"Final product contains {heterocycle}")

        # Return True if at least 2 different heterocycles are found
        return heterocycle_count >= 2
    except Exception as e:
        print(f"Error in heterocycle_rich_product_strategy: {e}")
        return False
