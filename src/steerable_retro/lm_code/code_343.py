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
    This function detects if the synthetic route involves a difluoromethyl group in the final product.
    """
    difluoromethyl_pattern = Chem.MolFromSmarts("[F][C]([F])[#6]")
    has_difluoromethyl = False

    def dfs_traverse(node, depth=0):
        nonlocal has_difluoromethyl

        if node["type"] == "mol" and depth == 0:  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(difluoromethyl_pattern):
                has_difluoromethyl = True
                print(f"Difluoromethyl group detected in final product: {node['smiles']}")

        # Traverse children with increased depth
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_difluoromethyl
