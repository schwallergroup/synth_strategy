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
    Detects a synthesis route that incorporates a 4-methoxyphenyl group
    via a sulfonamide formation.
    """
    has_methoxy_aryl_sulfonamide = False

    def dfs_traverse(node, depth=0):
        nonlocal has_methoxy_aryl_sulfonamide

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for 4-methoxyphenyl sulfonamide
                pattern = Chem.MolFromSmarts("[#7][S](=[O])(=[O])[c]1[c]([O][CH3])[c][c][c][c]1")
                if mol.HasSubstructMatch(pattern):
                    print("Detected 4-methoxyphenyl sulfonamide")
                    has_methoxy_aryl_sulfonamide = True

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    return has_methoxy_aryl_sulfonamide
