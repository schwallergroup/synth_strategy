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
    This function detects if the synthetic route involves the use of trimethylsilyl (TMS)
    protecting groups throughout the synthesis.
    """
    tms_found = False

    def dfs_traverse(node, depth=0):
        nonlocal tms_found

        if node["type"] == "mol":
            if "smiles" in node:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for TMS group
                    tms_pattern = Chem.MolFromSmarts("[C][Si]([C])([C])[C]")
                    if mol.HasSubstructMatch(tms_pattern):
                        print(f"Found TMS protecting group in molecule at depth {depth}")
                        tms_found = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return tms_found
