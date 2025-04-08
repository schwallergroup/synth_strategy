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
    This function detects if the synthesis retains a cyano group throughout multiple steps.
    """
    # Track paths with cyano groups
    cyano_paths = []

    def dfs_traverse(node, path=None, depth=0):
        if path is None:
            path = []

        current_path = path.copy()

        if node["type"] == "mol":
            mol_smiles = node["smiles"]
            has_cyano = checker.check_fg("Nitrile", mol_smiles)

            if has_cyano:
                print(f"Found cyano group at depth {depth}, SMILES: {mol_smiles}")
                current_path.append((depth, "mol", has_cyano))

                # If we have a path with at least 3 molecules containing cyano groups
                if len(current_path) >= 3:
                    # Check if these are connected through reactions (not just arbitrary molecules)
                    connected = True
                    for i in range(1, len(current_path)):
                        # Check if depths indicate connected molecules (should differ by 2 for mol->reaction->mol)
                        if current_path[i][0] - current_path[i - 1][0] != 2:
                            connected = False
                            break

                    if connected:
                        print(f"Found connected path with cyano groups: {current_path}")
                        cyano_paths.append(current_path)
                        return True
            else:
                # If this molecule doesn't have a cyano group, reset the path
                current_path = []

        # Continue traversal with children
        for child in node.get("children", []):
            if dfs_traverse(child, current_path, depth + 1):
                return True

        return False

    result = dfs_traverse(route)

    # If we found any valid paths during traversal
    if result or cyano_paths:
        return True

    return False
