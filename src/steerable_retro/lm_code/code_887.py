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
    Detects if the synthesis route is built around a thiophene core
    that remains present throughout the synthesis.
    """
    steps_with_thiophene = 0
    total_steps = 0

    def dfs_traverse(node, depth=0):
        nonlocal steps_with_thiophene, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check for thiophene ring
                thiophene_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#16][#6]1")

                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(thiophene_pattern):
                    steps_with_thiophene += 1
                    print(f"Thiophene core detected in product at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # Return True if thiophene is present in at least 75% of steps
    return total_steps > 0 and steps_with_thiophene / total_steps >= 0.75
