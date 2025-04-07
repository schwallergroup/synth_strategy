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
    This function detects if the final product contains a trifluoromethyl group.
    """
    has_trifluoromethyl = False

    def dfs_traverse(node):
        nonlocal has_trifluoromethyl

        if node["type"] == "mol" and node.get("smiles"):
            # Check if this is the final product (no parent reaction)
            if not any(
                child["type"] == "reaction" for child in node.get("children", [])
            ):
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    trifluoromethyl_pattern = Chem.MolFromSmarts("[C]([F])([F])[F]")
                    if mol.HasSubstructMatch(trifluoromethyl_pattern):
                        print("Detected trifluoromethyl group in final product")
                        has_trifluoromethyl = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return has_trifluoromethyl
