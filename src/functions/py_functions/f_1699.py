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
    This function detects a synthetic strategy involving a thiophene-containing scaffold.
    """
    has_thiophene_scaffold = False

    def dfs_traverse(node, depth=0):
        nonlocal has_thiophene_scaffold

        if node["type"] == "mol":
            smiles = node["smiles"]
            mol = Chem.MolFromSmiles(smiles)

            if mol:
                # Check for thiophene pattern
                thiophene_pattern = Chem.MolFromSmarts("c1cscc1")
                if mol.HasSubstructMatch(thiophene_pattern):
                    print(f"Found thiophene scaffold at depth {depth}")
                    has_thiophene_scaffold = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_thiophene_scaffold
