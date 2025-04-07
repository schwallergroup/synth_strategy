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
    Detects synthesis strategy involving protected nitrogen heterocycle.
    """
    protected_n_heterocycle = False

    def dfs_traverse(node):
        nonlocal protected_n_heterocycle

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for Boc-protected nitrogen in a ring
                boc_n_pattern = Chem.MolFromSmarts(
                    "[#6](=[O])-[O]-[C](-[CH3])(-[CH3])-[CH3]"
                )
                n_heterocycle_pattern = Chem.MolFromSmarts("[N;R]")

                if mol.HasSubstructMatch(boc_n_pattern) and mol.HasSubstructMatch(
                    n_heterocycle_pattern
                ):
                    print("Detected protected nitrogen heterocycle")
                    protected_n_heterocycle = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return protected_n_heterocycle
