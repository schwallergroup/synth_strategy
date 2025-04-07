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
    This function detects if the synthesis involves heterocyclic systems,
    specifically indole and thiazole/thiophene.
    """
    indole_present = False
    thiazole_or_thiophene_present = False

    def dfs_traverse(node, depth=0):
        nonlocal indole_present, thiazole_or_thiophene_present

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for indole
            if checker.check_ring("indole", mol_smiles):
                indole_present = True
                print(f"Indole detected in molecule: {mol_smiles}")

            # Check for thiazole or thiophene
            if checker.check_ring("thiazole", mol_smiles) or checker.check_ring(
                "thiophene", mol_smiles
            ):
                thiazole_or_thiophene_present = True
                print(f"Thiazole or thiophene detected in molecule: {mol_smiles}")

        # Check reaction nodes for thiazole formation
        elif node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check if the product contains thiazole or thiophene
            if checker.check_ring("thiazole", product) or checker.check_ring("thiophene", product):
                thiazole_or_thiophene_present = True
                print(f"Thiazole or thiophene detected in reaction product: {product}")

            # Check if the reaction is a thiazole formation
            if checker.check_reaction("benzothiazole", rsmi) or checker.check_reaction(
                "thiazole", rsmi
            ):
                thiazole_or_thiophene_present = True
                print(f"Thiazole formation reaction detected: {rsmi}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both heterocycles are present
    result = indole_present and thiazole_or_thiophene_present
    print(
        f"Final result: indole_present={indole_present}, thiazole_or_thiophene_present={thiazole_or_thiophene_present}, returning {result}"
    )
    return result
