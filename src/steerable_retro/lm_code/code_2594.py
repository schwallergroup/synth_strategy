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
    This function detects a synthetic strategy where nitrile functional groups
    are preserved throughout multiple steps of the synthesis.
    """
    nitrile_depths = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[C]#[N]")):
                    nitrile_depths.append(depth)
                    print(f"Detected nitrile at depth {depth}: {node['smiles']}")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Strategy is present if nitriles appear at multiple depths
    unique_depths = len(set(nitrile_depths))
    result = unique_depths >= 2
    print(
        f"Nitrile-preserving strategy detected: {result} (present at {unique_depths} different depths)"
    )
    return result
