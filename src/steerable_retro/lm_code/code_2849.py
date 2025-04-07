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
    This function detects if nitrile-containing intermediates are used in the synthesis.
    """
    nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
    nitrile_detected = False

    def dfs_traverse(node):
        nonlocal nitrile_detected

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(nitrile_pattern):
                print("Nitrile-containing intermediate detected")
                nitrile_detected = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return nitrile_detected
