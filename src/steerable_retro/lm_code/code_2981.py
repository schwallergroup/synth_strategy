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
    Detects if the synthesis route involves formation and use of hydrazide
    intermediates as key building blocks.
    """
    has_hydrazide_intermediate = False

    def dfs_traverse(node, depth=0):
        nonlocal has_hydrazide_intermediate

        if node["type"] == "mol" and depth > 0:
            # Check for hydrazide pattern in intermediates
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                hydrazide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][NX3]")
                if mol.HasSubstructMatch(hydrazide_pattern):
                    has_hydrazide_intermediate = True
                    print(f"Hydrazide intermediate found at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return has_hydrazide_intermediate
