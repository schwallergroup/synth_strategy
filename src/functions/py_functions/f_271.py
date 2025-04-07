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
    This function detects if the synthesis involves a proline derivative.
    """
    has_proline_scaffold = False

    # SMARTS pattern for proline scaffold
    proline_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#6][#6]1")

    def dfs_traverse(node):
        nonlocal has_proline_scaffold

        if node["type"] == "mol" and node.get("smiles"):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(proline_pattern):
                print("Found proline scaffold")
                has_proline_scaffold = True

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_proline_scaffold
