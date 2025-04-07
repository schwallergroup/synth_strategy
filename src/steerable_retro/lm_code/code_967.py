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
    This function detects if complex heterocyclic cores (pyrimidopyridazine) are preserved throughout the synthesis.
    """
    all_nodes_have_core = True

    def dfs_traverse(node, depth=0):
        nonlocal all_nodes_have_core

        if node["type"] == "mol" and not node.get("in_stock", False):
            # Check for pyrimidopyridazine core
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Simplified pattern for pyrimidopyridazine core
                core_pattern = Chem.MolFromSmarts("[n]1[c][c][n]2[c](=O)[c][c][n][c]12")
                if not mol.HasSubstructMatch(core_pattern):
                    all_nodes_have_core = False

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return all_nodes_have_core
