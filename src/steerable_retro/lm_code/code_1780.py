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
    This function detects a synthetic strategy involving a nitrofuran moiety
    in the final product.
    """
    has_nitrofuran = False

    def dfs_traverse(node, depth=0):
        nonlocal has_nitrofuran

        if node["type"] == "mol" and depth == 0:  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(
                Chem.MolFromSmarts("[#6]1[#6][#6][#6]([N+](=O)[O-])[#8]1")
            ):
                print("Detected nitrofuran moiety in final product")
                has_nitrofuran = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return has_nitrofuran
