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
    Detects if the synthesis route involves multiple different nitrogen-containing
    heterocycles (at least 3 different types).
    """
    heterocycle_patterns = {
        "benzothiazole": Chem.MolFromSmarts("c1ccc2c(c1)nc(s2)"),
        "thiazole": Chem.MolFromSmarts("c1nc(cs1)"),
        "pyrrolidine": Chem.MolFromSmarts("[#6]1[#6][#6][#7][#6]1"),
        "piperidine": Chem.MolFromSmarts("[#6]1[#6][#6][#7][#6][#6]1"),
        "tetrahydropyridine": Chem.MolFromSmarts("[#6]1[#6]=[#6][#7][#6][#6]1"),
    }

    found_heterocycles = set()

    def dfs_traverse(node):
        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                for name, pattern in heterocycle_patterns.items():
                    if mol.HasSubstructMatch(pattern):
                        found_heterocycles.add(name)

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)

    if len(found_heterocycles) >= 3:
        print(f"Found multiple nitrogen heterocycles: {found_heterocycles}")
        return True
    return False
