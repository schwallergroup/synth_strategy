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
    Detects if the synthesis involves a biaryl ether motif as a key structural element.
    """
    biaryl_ether_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal biaryl_ether_detected

        if node["type"] == "mol":
            # Check for biaryl ether pattern
            biaryl_ether_pattern = Chem.MolFromSmarts("[c][OX2][c]")
            mol = Chem.MolFromSmiles(node["smiles"])

            if mol and mol.HasSubstructMatch(biaryl_ether_pattern):
                print(f"Biaryl ether motif detected at depth {depth}")
                biaryl_ether_detected = True

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)
    return biaryl_ether_detected
