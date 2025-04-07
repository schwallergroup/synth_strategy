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
    Detects a synthesis involving a trifluoromethyl-containing compound.
    """
    has_trifluoromethyl = False

    def dfs_traverse(node):
        nonlocal has_trifluoromethyl

        if node["type"] == "mol" and node.get("smiles"):
            # Check for trifluoromethyl group
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                trifluoromethyl_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")
                if mol.HasSubstructMatch(trifluoromethyl_pattern):
                    has_trifluoromethyl = True
                    print("Detected trifluoromethyl group")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)

    return has_trifluoromethyl
