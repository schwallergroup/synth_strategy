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
    This function detects if the synthetic route involves fluorinated building blocks.
    """
    fluorinated_blocks = 0

    def dfs_traverse(node):
        nonlocal fluorinated_blocks

        if node["type"] == "mol" and node.get("in_stock", False):
            # Check starting materials for fluorine atoms
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                fluorine_pattern = Chem.MolFromSmarts("[F]")
                if mol.HasSubstructMatch(fluorine_pattern):
                    fluorinated_blocks += 1
                    print(f"Detected fluorinated building block: {node['smiles']}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    print(f"Total fluorinated building blocks detected: {fluorinated_blocks}")
    return fluorinated_blocks >= 1
