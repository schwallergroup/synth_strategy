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
    This function detects a strategy involving heterocyclic systems like piperazine
    and quinoline/pyridine derivatives.
    """
    piperazine_detected = False
    quinoline_detected = False

    def dfs_traverse(node):
        nonlocal piperazine_detected, quinoline_detected

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for piperazine
                piperazine_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6][#6]1")
                if mol.HasSubstructMatch(piperazine_pattern):
                    print("Found piperazine heterocycle")
                    piperazine_detected = True

                # Check for quinoline/pyridine systems
                quinoline_pattern = Chem.MolFromSmarts("[c]1[c][c][n][c]2[c][c][c][c][c]12")
                pyridine_pattern = Chem.MolFromSmarts("[c]1[c][c][n][c][c]1")

                if mol.HasSubstructMatch(quinoline_pattern) or mol.HasSubstructMatch(
                    pyridine_pattern
                ):
                    print("Found quinoline/pyridine heterocycle")
                    quinoline_detected = True

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    # Return True if both heterocycle types are detected
    return piperazine_detected and quinoline_detected
