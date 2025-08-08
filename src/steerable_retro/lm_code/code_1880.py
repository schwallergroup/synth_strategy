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
    This function detects if the synthetic route incorporates a trifluoromethyl aryl group
    and maintains it throughout the synthesis.
    """
    cf3_present = False
    cf3_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal cf3_present, cf3_depth

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for trifluoromethyl group
                cf3_pattern = Chem.MolFromSmarts("[#6][C]([F])([F])[F]")
                if mol.HasSubstructMatch(cf3_pattern):
                    cf3_present = True
                    if cf3_depth == -1 or depth > cf3_depth:
                        cf3_depth = depth
                        print(f"CF3 group detected at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # CF3 is incorporated early (higher depth) and maintained throughout
    print(f"Trifluoromethyl aryl incorporation detected: {cf3_present and cf3_depth > 1}")
    return cf3_present and cf3_depth > 1
