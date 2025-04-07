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
    Detects if a nitro-substituted aromatic fragment is used in the synthesis.
    """
    # Track nitro group presence
    has_nitro_aromatic = False

    # SMARTS pattern for nitro group on aromatic ring
    nitro_aromatic_pattern = "c-[N+](=[O])[O-]"

    def dfs_traverse(node, depth=0):
        nonlocal has_nitro_aromatic

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(Chem.MolFromSmarts(nitro_aromatic_pattern)):
                has_nitro_aromatic = True
                print(f"Detected nitro-aromatic group at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_nitro_aromatic:
        print("Detected nitro-aromatic containing fragment strategy")

    return has_nitro_aromatic
