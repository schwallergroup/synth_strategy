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
    Detects if the synthesis involves building or modifying a tetrahydronaphthalene core.
    """
    has_tetrahydronaphthalene = False

    def dfs_traverse(node):
        nonlocal has_tetrahydronaphthalene

        if node["type"] == "mol" and node.get("smiles"):
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # SMARTS pattern for tetrahydronaphthalene core
                    tetrahydronaphthalene_pattern = Chem.MolFromSmarts("c1cccc2c1CCCC2")
                    if mol.HasSubstructMatch(tetrahydronaphthalene_pattern):
                        has_tetrahydronaphthalene = True
                        print("Tetrahydronaphthalene core detected")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    return has_tetrahydronaphthalene
