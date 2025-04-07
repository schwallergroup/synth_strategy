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
    This function detects a linear synthesis strategy focused on heterocycle modification.
    """
    # Track modifications to heterocycle
    modifications = []
    has_heterocycle = False

    def dfs_traverse(node):
        nonlocal modifications, has_heterocycle

        # Check for heterocycle presence
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for heterocyclic structures
                heterocycle_pattern = Chem.MolFromSmarts("c1[nH]cc*1")  # Simplified pattern
                if mol.HasSubstructMatch(heterocycle_pattern):
                    has_heterocycle = True

        # Check for reactions modifying the heterocycle
        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            depth = node["metadata"].get("depth", 0)

            # Track the modification with its depth
            modifications.append((rsmi, depth))

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Sort modifications by depth to check if they're sequential
    modifications.sort(key=lambda x: x[1], reverse=True)

    # Check if we have at least 3 modifications and they're sequential
    is_linear = len(modifications) >= 3 and has_heterocycle

    if is_linear:
        print(f"Detected linear heterocycle modification with {len(modifications)} steps")

    return is_linear
