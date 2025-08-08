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
    This function detects if the route maintains a difluoro group throughout the synthesis.
    """
    difluoro_present_in_final = False
    difluoro_present_in_intermediates = False

    def dfs_traverse(node):
        nonlocal difluoro_present_in_final, difluoro_present_in_intermediates

        if node["type"] == "mol":
            # Check for difluoro group
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                difluoro_pattern = Chem.MolFromSmarts("[#6]([F])([F])")
                if mol.HasSubstructMatch(difluoro_pattern):
                    if node.get("in_stock", False):
                        # This is a starting material
                        difluoro_present_in_intermediates = True
                    elif not node.get("children"):
                        # This is the final product
                        difluoro_present_in_final = True
                    else:
                        # This is an intermediate
                        difluoro_present_in_intermediates = True

        # Continue traversing
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # The strategy is present if difluoro is in both final product and intermediates
    if difluoro_present_in_final and difluoro_present_in_intermediates:
        print("Detected maintenance of difluoro group throughout synthesis")
        return True
    return False
