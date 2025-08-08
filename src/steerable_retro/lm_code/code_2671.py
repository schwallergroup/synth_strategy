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
    Detects if the synthesis route involves a trifluoromethyl-containing aromatic fragment
    that persists throughout the synthesis.
    """
    depths_with_trifluoromethyl = set()
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal max_depth
        max_depth = max(max_depth, depth)

        if node["type"] == "mol" and node.get("smiles"):
            # Check for trifluoromethyl aromatic pattern
            mol = Chem.MolFromSmiles(node["smiles"])
            trifluoromethyl_pattern = Chem.MolFromSmarts("c-C(F)(F)F")

            if mol and trifluoromethyl_pattern and mol.HasSubstructMatch(trifluoromethyl_pattern):
                depths_with_trifluoromethyl.add(depth)
                print(f"Trifluoromethyl aromatic detected at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if trifluoromethyl is present at multiple depths
    result = len(depths_with_trifluoromethyl) >= 3 and 0 in depths_with_trifluoromethyl
    print(f"Trifluoromethyl aromatic persistence strategy detected: {result}")
    return result
