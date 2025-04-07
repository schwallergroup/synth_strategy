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
    Detects if a trifluoromethyl group on an aromatic ring is maintained
    throughout the synthesis.
    """
    cf3_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for CF3 group on aromatic ring
                    cf3_pattern = Chem.MolFromSmarts("c[C]([F])([F])[F]")
                    if mol.HasSubstructMatch(cf3_pattern):
                        cf3_depths.append(depth)
                        print(f"Found CF3 group at depth {depth}")
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if CF3 group is present at multiple depths (maintained throughout synthesis)
    return len(cf3_depths) >= 2
