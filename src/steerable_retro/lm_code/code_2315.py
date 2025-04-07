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
    This function detects if a piperidine scaffold is maintained throughout the synthesis.
    """
    piperidine_in_final = False
    piperidine_in_intermediates = []

    def dfs_traverse(node, depth=0):
        nonlocal piperidine_in_final, piperidine_in_intermediates

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            piperidine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#6][#6][#6]1")

            if mol and mol.HasSubstructMatch(piperidine_pattern):
                if depth == 0:
                    piperidine_in_final = True
                else:
                    piperidine_in_intermediates.append(depth)

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if piperidine is present in final product and in at least one intermediate
    if piperidine_in_final and piperidine_in_intermediates:
        print(f"Detected piperidine scaffold maintenance throughout synthesis")
        return True
    return False
