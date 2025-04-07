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
    Detects if the 1,3-benzodioxole (methylenedioxy) scaffold is preserved throughout the synthesis.
    """
    methylenedioxy_pattern = Chem.MolFromSmarts("c1cc2OCOc2cc1")
    all_steps_have_pattern = True

    def dfs_traverse(node):
        nonlocal all_steps_have_pattern

        if node["type"] == "mol" and node["smiles"]:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol is not None and not mol.HasSubstructMatch(methylenedioxy_pattern):
                if not node.get("in_stock", False):  # Ignore starting materials
                    all_steps_have_pattern = False
                    print(f"Molecule without methylenedioxy pattern: {node['smiles']}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return all_steps_have_pattern
