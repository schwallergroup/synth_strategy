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
    Detects if the synthesis route uses an alkyne as a linking group
    between aromatic systems.
    """
    alkyne_linker_detected = False

    def dfs_traverse(node):
        nonlocal alkyne_linker_detected

        if node["type"] == "mol":
            smiles = node.get("smiles", "")
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Check for alkyne connecting two aromatic systems
                    alkyne_linker_pattern = Chem.MolFromSmarts("c[C]#[C]c")
                    if mol.HasSubstructMatch(alkyne_linker_pattern):
                        print("Alkyne linker detected")
                        alkyne_linker_detected = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return alkyne_linker_detected
