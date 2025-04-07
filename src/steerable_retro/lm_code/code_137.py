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
    This function detects if the synthesis involves heterocyclic compounds,
    particularly focusing on pyrazole and other nitrogen-containing heterocycles.
    """
    has_heterocycle = False

    def dfs_traverse(node):
        nonlocal has_heterocycle

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for pyrazole pattern
                pyrazole_pattern = Chem.MolFromSmarts("n1ncc[c]1")
                # Check for other nitrogen heterocycles
                n_heterocycle_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#6][#6]1")

                if mol.HasSubstructMatch(pyrazole_pattern) or mol.HasSubstructMatch(
                    n_heterocycle_pattern
                ):
                    print(f"Found heterocycle in molecule: {node['smiles']}")
                    has_heterocycle = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_heterocycle
