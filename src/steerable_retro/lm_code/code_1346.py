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
    Detects if the synthesis route includes BOC-protected amines.
    """
    boc_protection_found = False

    def dfs_traverse(node):
        nonlocal boc_protection_found

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                boc_pattern = Chem.MolFromSmarts("[CH3]C([CH3])([CH3])[O]C(=O)[NH]")
                if mol.HasSubstructMatch(boc_pattern):
                    boc_protection_found = True
                    print("BOC protecting group detected")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return boc_protection_found
