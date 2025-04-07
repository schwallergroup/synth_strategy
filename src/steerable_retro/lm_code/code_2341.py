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
    throughout the synthesis.
    """
    steps_with_diaryl_ether = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_diaryl_ether, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for diaryl ether pattern
                diaryl_ether_pattern = Chem.MolFromSmarts("[c]-[O]-[c]")

                product_mol = Chem.MolFromSmiles(product) if product else None

                if product_mol and product_mol.HasSubstructMatch(diaryl_ether_pattern):
                    print("Found diaryl ether at depth", depth)
                    steps_with_diaryl_ether += 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    # Return True if diaryl ether is present in all steps
    return steps_with_diaryl_ether > 0 and steps_with_diaryl_ether == total_steps
