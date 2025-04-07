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
    Detects a synthesis route that includes formation of an isothiocyanate intermediate.
    """
    has_isothiocyanate_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_isothiocyanate_formation

        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            product_part = rsmi.split(">")[-1]

            # Check for isothiocyanate formation
            isothiocyanate_pattern = Chem.MolFromSmarts("[#7]=[#6]=[#16]")
            product_mol = Chem.MolFromSmiles(product_part)

            if product_mol and product_mol.HasSubstructMatch(isothiocyanate_pattern):
                print(f"Found isothiocyanate formation at depth {depth}")
                has_isothiocyanate_formation = True

        # Continue traversing the tree
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    if has_isothiocyanate_formation:
        print("Isothiocyanate formation sequence detected")
    else:
        print("Isothiocyanate formation sequence NOT detected")

    return has_isothiocyanate_formation
