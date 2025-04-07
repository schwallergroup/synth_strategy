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
    This function detects if the synthesis involves quinoline-containing
    compounds as key intermediates or products.
    """
    quinoline_found = False

    def dfs_traverse(node, depth=0):
        nonlocal quinoline_found

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol_smiles = node["smiles"]
                mol = Chem.MolFromSmiles(mol_smiles)

                if mol:
                    # Check for quinoline structure using the checker function
                    if checker.check_ring("quinoline", mol_smiles):
                        print(f"Quinoline ring detected at depth {depth}, SMILES: {mol_smiles}")
                        quinoline_found = True

                    # Also check for isoquinoline which is a structural isomer
                    elif checker.check_ring("isoquinoline", mol_smiles):
                        print(f"Isoquinoline ring detected at depth {depth}, SMILES: {mol_smiles}")
                        quinoline_found = True
                else:
                    print(
                        f"Warning: Could not parse molecule SMILES at depth {depth}: {mol_smiles}"
                    )
            except Exception as e:
                print(f"Error processing molecule at depth {depth}: {e}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Quinoline containing strategy result: {quinoline_found}")
    return quinoline_found
