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
    This function detects if the synthetic route involves fluorinated aromatic building blocks.
    """
    has_fluorinated_aromatics = False

    def dfs_traverse(node):
        nonlocal has_fluorinated_aromatics

        if node["type"] == "mol" and node.get("in_stock", False):
            # Check starting materials for fluorinated aromatics
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                fluoro_aromatic_pattern = Chem.MolFromSmarts("c[F]")
                if mol.HasSubstructMatch(fluoro_aromatic_pattern):
                    has_fluorinated_aromatics = True
                    print(
                        f"Found fluorinated aromatic building block: {node['smiles']}"
                    )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Fluorinated aromatics strategy detected: {has_fluorinated_aromatics}")
    return has_fluorinated_aromatics
