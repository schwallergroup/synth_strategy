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
    Detects synthesis involving a fluorophenyl group.
    """
    fluorophenyl_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal fluorophenyl_detected

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Fluorophenyl pattern
                fluorophenyl_pattern = Chem.MolFromSmarts("c1ccccc1F")
                if mol.HasSubstructMatch(fluorophenyl_pattern):
                    fluorophenyl_detected = True
                    print(f"Detected fluorophenyl group at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return fluorophenyl_detected
