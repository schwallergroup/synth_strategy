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
    Detects if the synthesis route maintains a diaryl ether linkage throughout.
    """
    steps_with_diaryl_ether = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal steps_with_diaryl_ether, total_steps

        if node["type"] == "reaction":
            total_steps += 1

            # Get reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            product = rsmi.split(">")[-1]

            # Check for diaryl ether pattern in product
            diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")
            product_mol = Chem.MolFromSmiles(product)

            if product_mol is not None and product_mol.HasSubstructMatch(diaryl_ether_pattern):
                steps_with_diaryl_ether += 1
                print(f"Detected diaryl ether in step {total_steps}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if diaryl ether is present in all or most steps (>80%)
    return total_steps > 0 and steps_with_diaryl_ether / total_steps > 0.8
