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
    This function detects if a heterocyclic aromatic core (thiazole) is preserved throughout the synthesis.
    """
    # Track if thiazole is present in all molecules
    all_mols_have_thiazole = True

    def dfs_traverse(node):
        nonlocal all_mols_have_thiazole

        if node["type"] == "mol" and "smiles" in node:
            # Skip starting materials
            if not node.get("in_stock", False):
                # Check if this molecule has a thiazole ring using the checker function
                has_thiazole = checker.check_ring("thiazole", node["smiles"])

                if not has_thiazole:
                    # If any molecule doesn't have thiazole, set flag to False
                    all_mols_have_thiazole = False
                    print(f"Molecule without thiazole found: {node['smiles']}")
                else:
                    print(f"Molecule with thiazole found: {node['smiles']}")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return all_mols_have_thiazole
