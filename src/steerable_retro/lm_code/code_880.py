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
    This function detects a synthetic strategy involving a gem-difluoro motif
    that is maintained throughout the synthesis.
    """
    gem_difluoro_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal gem_difluoro_count

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None:
                gem_difluoro_pattern = Chem.MolFromSmarts("[CX4]([FX1])([FX1])")
                if mol.HasSubstructMatch(gem_difluoro_pattern):
                    gem_difluoro_count += 1
                    print(f"Found gem-difluoro motif at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # If gem-difluoro appears in multiple steps, it's maintained throughout
    return gem_difluoro_count >= 2
