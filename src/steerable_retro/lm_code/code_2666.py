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
    Detects if the synthesis involves a cyclopropyl group, particularly as part of a sulfonamide.
    """
    has_cyclopropyl = False

    def dfs_traverse(node):
        nonlocal has_cyclopropyl

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for cyclopropyl group
                cyclopropyl_pattern = Chem.MolFromSmarts("C1CC1")
                if mol.HasSubstructMatch(cyclopropyl_pattern):
                    # Check if it's part of a sulfonamide
                    cyclopropyl_sulfonamide = Chem.MolFromSmarts("[NH]S(=O)(=O)C1CC1")
                    if mol.HasSubstructMatch(cyclopropyl_sulfonamide):
                        print("Found cyclopropyl sulfonamide in molecule")
                        has_cyclopropyl = True
                    else:
                        print("Found cyclopropyl group in molecule")
                        has_cyclopropyl = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return has_cyclopropyl
