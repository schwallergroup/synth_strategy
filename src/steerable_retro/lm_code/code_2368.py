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
    This function detects if the synthesis maintains halogenated aromatics throughout.
    """
    chloro_aromatic_pattern = Chem.MolFromSmarts("c[Cl]")
    depths_with_chloro = set()

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(chloro_aromatic_pattern):
                depths_with_chloro.add(depth)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if chlorinated aromatics are present at multiple depths
    if len(depths_with_chloro) >= 3:
        print(
            f"Halogenated aromatics maintained throughout synthesis at depths: {depths_with_chloro}"
        )
        return True
    return False
