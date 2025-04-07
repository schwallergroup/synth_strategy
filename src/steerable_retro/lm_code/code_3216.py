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
    Detects if the synthesis route involves compounds containing
    methoxy-substituted aromatic rings.
    """
    has_methoxy_aromatic = False

    def dfs_traverse(node, depth=0):
        nonlocal has_methoxy_aromatic

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for methoxy aromatic
                methoxy_aromatic_pattern = Chem.MolFromSmarts("[c]-[O][C]")
                if mol.HasSubstructMatch(methoxy_aromatic_pattern):
                    print(f"Found methoxy-substituted aromatic ring at depth {depth}")
                    has_methoxy_aromatic = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)
    return has_methoxy_aromatic
