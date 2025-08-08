#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter


def main(route):
    """
    Detects a strategy involving the use of a diphenylmethyl fragment throughout the synthesis.
    """
    diphenylmethyl_present = False

    def dfs_traverse(node):
        nonlocal diphenylmethyl_present

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # SMARTS pattern for diphenylmethyl group
                diphenylmethyl_pattern = Chem.MolFromSmarts(
                    "[#6](-[c]1[c][c][c][c][c]1)(-[c]1[c][c][c][c][c]1)"
                )
                if mol.HasSubstructMatch(diphenylmethyl_pattern):
                    diphenylmethyl_present = True
                    print("Found diphenylmethyl group in molecule")

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    print(f"Diphenylmethyl fragment strategy detected: {diphenylmethyl_present}")
    return diphenylmethyl_present
