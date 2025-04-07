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
    Detects if the synthetic route involves a halogenated heterocycle (specifically fluoroindazole).
    """
    contains_fluoroindazole = False

    def dfs_traverse(node):
        nonlocal contains_fluoroindazole

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for fluoroindazole pattern
                fluoroindazole_pattern = Chem.MolFromSmarts("c1cc2[nH]nc(F)c2cc1")
                if mol.HasSubstructMatch(fluoroindazole_pattern):
                    print("Fluoroindazole detected")
                    contains_fluoroindazole = True

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return contains_fluoroindazole
