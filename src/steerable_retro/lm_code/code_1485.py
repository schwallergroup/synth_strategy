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
    This function detects a strategy involving the incorporation of a
    trifluoromethyl group in the synthesis.
    """
    has_trifluoromethyl = False

    def dfs_traverse(node):
        nonlocal has_trifluoromethyl

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for trifluoromethyl group
                    cf3_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")
                    if mol.HasSubstructMatch(cf3_pattern):
                        has_trifluoromethyl = True
                        print(f"Found trifluoromethyl group in molecule: {node['smiles']}")
            except Exception as e:
                print(f"Error in trifluoromethyl detection: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Trifluoromethyl group present: {has_trifluoromethyl}")
    return has_trifluoromethyl
