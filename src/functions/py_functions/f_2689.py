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
    This function detects preservation of indole heterocycle throughout the synthesis.
    """
    indole_present_all_steps = True
    indole_pattern = Chem.MolFromSmarts(
        "[#6]1:[#6]:[#6]:[#6]2:[#7]:[#6]:[#6]:[#6]:[#6]:2:[#6]:1"
    )

    def dfs_traverse(node, depth=0):
        nonlocal indole_present_all_steps

        if node["type"] == "mol" and not node.get("in_stock", False):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and not mol.HasSubstructMatch(indole_pattern):
                indole_present_all_steps = False

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return indole_present_all_steps
