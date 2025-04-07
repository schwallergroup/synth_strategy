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
    Detects if the synthesis involves a methoxy-substituted aromatic ring as a key structural element.
    """
    has_methoxy_aryl = False

    def dfs_traverse(node):
        nonlocal has_methoxy_aryl

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for methoxy on aromatic ring
                methoxy_aryl_pattern = Chem.MolFromSmarts("c[OX2][CH3]")
                if mol.HasSubstructMatch(methoxy_aryl_pattern):
                    print("Found methoxy-substituted aromatic ring")
                    has_methoxy_aryl = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from root
    dfs_traverse(route)
    return has_methoxy_aryl
