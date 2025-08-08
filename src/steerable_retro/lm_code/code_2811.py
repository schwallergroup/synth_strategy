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
    Detects if the synthesis involves an indole-containing building block.
    """
    indole_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal indole_detected

        if node["type"] == "mol":
            smiles = node["smiles"]
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # Indole pattern
                    indole_pattern = Chem.MolFromSmarts("[#6]1[#6][#7H][#6]2[#6][#6][#6][#6][#6]12")
                    if mol.HasSubstructMatch(indole_pattern):
                        print(f"Detected indole-containing structure at depth {depth}")
                        indole_detected = True
            except:
                print("Error processing molecule SMILES")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)
    return indole_detected
