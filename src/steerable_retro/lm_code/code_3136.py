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
    Detects synthesis routes that maintain halogen functionality (particularly iodine)
    throughout the synthesis, suggesting it's being preserved for later functionalization.
    """
    # Track halogen presence at each step
    halogen_at_depths = {}
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node.get("type") == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Check for iodine
                    iodine_pattern = Chem.MolFromSmarts("[#53]")
                    # Check for any halogen
                    halogen_pattern = Chem.MolFromSmarts("[#9,#17,#35,#53]")

                    if mol.HasSubstructMatch(iodine_pattern):
                        halogen_at_depths[depth] = "I"
                    elif mol.HasSubstructMatch(halogen_pattern):
                        halogen_at_depths[depth] = "X"  # Other halogen
            except:
                pass

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if halogen (particularly iodine) is preserved throughout
    if len(halogen_at_depths) < 2:
        return False

    # Check if iodine is present at the beginning (highest depth) and end (depth 0)
    has_iodine_at_start = halogen_at_depths.get(max_depth) == "I"
    has_iodine_at_end = halogen_at_depths.get(0) == "I"

    if has_iodine_at_start and has_iodine_at_end:
        print("Found halogen (iodine) preservation throughout synthesis")
        return True

    return False
