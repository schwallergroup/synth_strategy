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
    This function detects if the synthesis involves fluorinated substituents throughout.
    """
    fluorine_count = 0
    total_steps = 0

    def dfs_traverse(node):
        nonlocal fluorine_count, total_steps

        if node["type"] == "reaction":
            total_steps += 1
            if "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                product = rsmi.split(">")[-1]

                # Check if product contains fluorine
                product_mol = Chem.MolFromSmiles(product)
                if product_mol and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[F]")
                ):
                    fluorine_count += 1
                    print("Found fluorinated intermediate")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # If most steps (>80%) involve fluorinated compounds, return True
    return total_steps > 0 and (fluorine_count / total_steps) > 0.8
