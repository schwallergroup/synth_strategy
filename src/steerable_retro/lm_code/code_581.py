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
    Detects if a trifluoromethyl group is introduced in the late stages of synthesis.
    """
    trifluoromethyl_depth = None
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal trifluoromethyl_depth, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("C(F)(F)F")):
                if trifluoromethyl_depth is None or depth < trifluoromethyl_depth:
                    trifluoromethyl_depth = depth

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if CF3 is introduced in the late stage (first half of synthesis)
    # Remember that low depth = late stage in synthesis
    is_late_stage = trifluoromethyl_depth is not None and trifluoromethyl_depth < (max_depth / 2)

    print(f"Late stage trifluoromethyl introduction: {is_late_stage}")
    print(f"Trifluoromethyl first appears at depth: {trifluoromethyl_depth}")
    print(f"Maximum synthesis depth: {max_depth}")

    return is_late_stage
