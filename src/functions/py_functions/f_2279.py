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
    Detects if the synthesis involves pyrimidine-based intermediates.
    """
    found_pyrimidine = False

    def dfs_traverse(node):
        nonlocal found_pyrimidine

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Pattern for pyrimidine
                pyrimidine_pattern = Chem.MolFromSmarts(
                    "[#6]1:[#7]:[#6]:[#6]:[#7]:[#6]1"
                )

                if mol.HasSubstructMatch(pyrimidine_pattern):
                    print("Found pyrimidine-containing intermediate")
                    found_pyrimidine = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)
    return found_pyrimidine
