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


def main(route):
    """
    This function detects if an indole core is maintained throughout the synthesis.
    """
    indole_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#6]2[#6][#6][#6][#6][#6]12")
    indole_present_at_all_depths = True

    def dfs_traverse(node):
        nonlocal indole_present_at_all_depths
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                if not mol.HasSubstructMatch(indole_pattern):
                    if not node.get(
                        "in_stock", False
                    ):  # Only check non-starting materials
                        indole_present_at_all_depths = False
                        print(f"Indole core not found in molecule: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Indole core maintained throughout: {indole_present_at_all_depths}")
    return indole_present_at_all_depths
