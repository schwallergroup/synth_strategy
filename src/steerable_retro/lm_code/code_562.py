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
    This function detects if the synthesis involves a piperidine-containing compound
    throughout the synthesis.
    """
    piperidine_present = False

    def dfs_traverse(node, depth=0):
        nonlocal piperidine_present

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Piperidine pattern
                piperidine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#6][#6][#6]1")
                if mol.HasSubstructMatch(piperidine_pattern):
                    print(f"Found piperidine at depth {depth}")
                    piperidine_present = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return piperidine_present
