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
    This function detects if the synthetic route incorporates pyridine-containing building blocks.
    """
    has_pyridine = False

    def dfs_traverse(node):
        nonlocal has_pyridine

        if node["type"] == "mol" and node.get("in_stock", False):
            # Check starting materials for pyridine structures
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")
                if mol.HasSubstructMatch(pyridine_pattern):
                    has_pyridine = True
                    print(f"Found pyridine-containing building block: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Pyridine building blocks detected: {has_pyridine}")
    return has_pyridine
