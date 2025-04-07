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
    This function detects if the synthesis maintains a diaryl ether linkage
    throughout the synthetic route.
    """
    steps_with_diaryl_ether = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal steps_with_diaryl_ether, total_steps

        if node["type"] == "reaction" and node.get("metadata", {}).get("rsmi"):
            total_steps += 1
            rsmi = node["metadata"]["rsmi"]
            product_part = rsmi.split(">")[-1]

            # Pattern for diaryl ether
            diaryl_ether_pattern = Chem.MolFromSmarts("[c][O][c]")

            try:
                product_mol = Chem.MolFromSmiles(product_part)
                if product_mol and product_mol.HasSubstructMatch(diaryl_ether_pattern):
                    steps_with_diaryl_ether += 1
                    print(f"Detected diaryl ether in step {total_steps}")
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if diaryl ether is maintained throughout (in at least 80% of steps)
    return total_steps > 0 and (steps_with_diaryl_ether / total_steps) >= 0.8
