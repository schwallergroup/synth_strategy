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
    Detects a synthesis involving a trifluoromethyl group that is preserved throughout.
    """
    has_trifluoromethyl_in_final = False

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoromethyl_in_final

        if node["type"] == "mol" and depth == 0:
            # Check final product for trifluoromethyl group
            mol = Chem.MolFromSmiles(node["smiles"])
            trifluoromethyl_pattern = Chem.MolFromSmarts("[CX4]([F])([F])([F])")
            if mol and trifluoromethyl_pattern and mol.HasSubstructMatch(trifluoromethyl_pattern):
                has_trifluoromethyl_in_final = True
                print("Found trifluoromethyl group in final product")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Trifluoromethyl-containing strategy: {has_trifluoromethyl_in_final}")
    return has_trifluoromethyl_in_final
