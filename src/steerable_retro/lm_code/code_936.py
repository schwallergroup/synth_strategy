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
    This function detects if the synthesis involves triflate as an activating group.
    """
    triflate_pattern = Chem.MolFromSmarts("[#8]S(=O)(=O)[#6]([F])([F])[F]")

    triflate_used = False

    def dfs_traverse(node):
        nonlocal triflate_used

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(triflate_pattern):
                triflate_used = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Triflate activation strategy: {triflate_used}")
    return triflate_used
