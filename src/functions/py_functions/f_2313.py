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
    This function detects if a nitrile group is preserved throughout the synthesis.
    """
    nitrile_in_final = False
    nitrile_in_intermediates = []

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_in_final, nitrile_in_intermediates

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")

            if mol and mol.HasSubstructMatch(nitrile_pattern):
                if depth == 0:
                    nitrile_in_final = True
                else:
                    nitrile_in_intermediates.append(depth)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if nitrile is present in final product and in at least one intermediate
    if nitrile_in_final and nitrile_in_intermediates:
        print(f"Detected nitrile preservation throughout synthesis")
        return True
    return False
