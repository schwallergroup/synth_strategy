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
    Detects if the synthetic route involves a pyrazole-containing fragment
    """
    found_pyrazole = False

    def dfs_traverse(node):
        nonlocal found_pyrazole

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                # SMARTS pattern for pyrazole
                pyrazole_pattern = Chem.MolFromSmarts("c1nn([#6])cc1")
                if mol and mol.HasSubstructMatch(pyrazole_pattern):
                    print(f"Found pyrazole-containing fragment: {node['smiles']}")
                    found_pyrazole = True
            except:
                pass

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return found_pyrazole
