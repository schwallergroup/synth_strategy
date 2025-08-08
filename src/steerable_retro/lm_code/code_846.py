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
    This function detects a synthetic strategy involving a cyclohexyl group maintained
    throughout the synthesis.
    """
    # Track if cyclohexyl is present in final product
    cyclohexyl_in_final = False
    cyclohexyl_in_intermediates = 0

    def dfs_traverse(node, depth=0):
        nonlocal cyclohexyl_in_final, cyclohexyl_in_intermediates

        if node["type"] == "mol":
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    cyclohexyl_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]1")
                    if mol.HasSubstructMatch(cyclohexyl_pattern):
                        if depth == 0:  # Final product
                            cyclohexyl_in_final = True
                            print("Found cyclohexyl in final product")
                        else:
                            cyclohexyl_in_intermediates += 1
                            print(f"Found cyclohexyl in intermediate at depth {depth}")
            except:
                pass

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return cyclohexyl_in_final and cyclohexyl_in_intermediates > 0
