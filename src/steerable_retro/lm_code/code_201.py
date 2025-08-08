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
    This function detects if a pyridine motif is present throughout the synthesis.
    """
    # Track if pyridine is present at each depth
    pyridine_at_depth = {}

    def dfs_traverse(node, depth=0):
        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    pyridine_pattern = Chem.MolFromSmarts("[n]1[c][c][c][c][c]1")
                    if mol.HasSubstructMatch(pyridine_pattern):
                        pyridine_at_depth[depth] = True
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if pyridine is present at multiple depths
    result = len(pyridine_at_depth) >= 2
    if result:
        print(f"Pyridine motif found at multiple depths: {list(pyridine_at_depth.keys())}")

    return result
