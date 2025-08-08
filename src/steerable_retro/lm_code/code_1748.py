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
    This function detects if the synthesis route involves a tetralin scaffold
    that is maintained throughout the synthesis.
    """
    tetralin_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal tetralin_count

        if node["type"] == "mol" and "smiles" in node:
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Tetralin scaffold pattern
                tetralin_pattern = Chem.MolFromSmarts(
                    "[#6]1[#6][#6][#6]2[#6][#6][#6][#6][#6]2[#6]1"
                )
                if mol.HasSubstructMatch(tetralin_pattern):
                    tetralin_count += 1
                    print(f"Detected tetralin scaffold at depth {depth}")

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    # Return True if tetralin scaffold appears in multiple intermediates
    return tetralin_count >= 2
